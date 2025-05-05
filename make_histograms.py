from array import array
import os
import numpy as np
import ROOT as pyr
import yaml
import argparse

PYROOT_DEFAULT_DIR = pyr.gDirectory.pwd()

parser = argparse.ArgumentParser()
parser.add_argument("yamlpath", help="YAML spec file path")
parser.add_argument("--diagnosis", help="Create diagnosis file, showing event contributions from each input ROOT file", action="store_true")
parser.add_argument("--debug", help="Prints debugging messages", action="store_true")
args = parser.parse_args()

def print_histogram(hist):
    for i in range(1, hist.GetNbinsX()+1): print(f"{hist.GetBinContent(i):.2f}", end=' ')
    print()

class AnalysisHistogram(object):
    def __init__(self, categories, xbins, xmin, xmax, treename):
        self.categories = categories
        self.xbins = xbins
        self.xmin = xmin
        self.xmax = xmax
        self.nom_hist = {}
        self.unc_hist = {}
        self.data_hist = {}
        for category in self.categories:
            self.nom_hist[category] = {}
            self.unc_hist[category] = {}
    
    def add_nominal_hist(self, category, histobj, isPass=True):
        if category not in self.nom_hist.keys():
            raise KeyError(f"{category} is not defined. The defined categories are {self.categories}")
        self.nom_hist[category]["pass" if isPass else "fail"] = histobj
    
    def define_unc(self, category, unc):
        self.unc_hist[category][unc] = {}
        self.unc_hist[category][unc]["up"] = {"pass": {}, "fail": {}}
        self.unc_hist[category][unc]["down"] = {"pass": {}, "fail": {}}
    
    def add_unc_hist(self, category, unc, histobj, isUp=True, isPass=True):
        self.unc_hist[category][unc]["up" if isUp else "down"]["pass" if isPass else "fail"] = histobj.Clone()
    
    def add_data_hist(self, hist, isPass=True):
        self.data_hist["pass" if isPass else "fail"] = hist
    
    def save_histograms(self, filename):
        self.check_zero_bins()
        savefile = pyr.TFile(filename, "RECREATE")
        for category in self.categories:
            self.nom_hist[category]["pass"].Write()
            self.nom_hist[category]["fail"].Write()
            for unc in self.unc_hist[category].keys():
                self.unc_hist[category][unc]["up"]["pass"].Write()
                self.unc_hist[category][unc]["up"]["fail"].Write()
                self.unc_hist[category][unc]["down"]["pass"].Write()
                self.unc_hist[category][unc]["down"]["fail"].Write()
        self.data_hist["pass"].Write()
        self.data_hist["fail"].Write()
        savefile.Close()
        
    def check_zero_bins(self):
        for category in self.categories:
            for passing in ["pass", "fail"]:
                for i in range(1, self.nom_hist[category][passing].GetNbinsX()+1): 
                    if self.nom_hist[category][passing].GetBinContent(i) <=0:
                        self.nom_hist[category][passing].SetBinContent(i, 0.01)
                        self.nom_hist[category][passing].SetBinError(i, 0.01)
                for unc in self.unc_hist[category].keys():
                    for unctype in ["up", "down"]:
                        for i in range(1, self.unc_hist[category][unc][unctype][passing].GetNbinsX()+1): 
                            if self.unc_hist[category][unc][unctype][passing].GetBinContent(i) <=0:
                                self.unc_hist[category][unc][unctype][passing].SetBinContent(i, 0.01)
                                self.unc_hist[category][unc][unctype][passing].SetBinError(i, 0.01)

def extract_histogram(filename, treename, var, cut, weight, histname, xbins, xmin, xmax):
    fileobj = pyr.TFile(filename, "READ")
    treeobj = fileobj.Get(treename)
    hist = pyr.TH1F(histname, histname, xbins, xmin, xmax)
    if args.debug: print(f"({cut})*({weight})")
    project_out = treeobj.Project(histname, var, f"({cut})*({weight})", "e")
    if args.debug:
        print(project_out)
        print_histogram(hist)
    hist.SetDirectory(pyr.gROOT)
    fileobj.Close()
    return hist

def combine_histograms(histlist, finalname, xbins, xmin, xmax):
    cachehist = pyr.TH1F(finalname, finalname, xbins, xmin, xmax)
    if args.debug: print(histlist)
    for hist in histlist: cachehist.Add(hist)
    return cachehist

def create_lowess(nominal, upper, lower, symmetrize=True, limit=0.5):
    x_edges = [nominal.GetBinLowEdge(i) for i in range(1, nominal.GetNbinsX()+2)]
    x_points = [(x_edges[i]+x_edges[i+1])/2 for i in range(len(x_edges)-1)]
    
    nominal_y = np.array([nominal.GetBinContent(i) for i in range(1, nominal.GetNbinsX()+1)])
    upper_y   = np.array([upper.GetBinContent(i) for i in range(1, upper.GetNbinsX()+1)])
    lower_y   = np.array([lower.GetBinContent(i) for i in range(1, lower.GetNbinsX()+1)])

    graph_upper = pyr.TGraphErrors()
    graph_lower = pyr.TGraphErrors()
    for i, index in enumerate(np.flatnonzero(nominal_y)):
        graph_upper.SetPoint(i, x_points[index], 1+(upper_y[index]-lower_y[index])/nominal_y[index] if symmetrize else upper_y[index]/nominal_y[index])
        graph_lower.SetPoint(i, x_points[index], 1-(upper_y[index]-lower_y[index])/nominal_y[index] if symmetrize else lower_y[index]/nominal_y[index])

    graph_upper_smooth = pyr.TGraphSmooth()
    graph_lower_smooth = pyr.TGraphSmooth()

    graph_upper_out = graph_upper_smooth.SmoothLowess(graph_upper, "", 1.0)
    graph_lower_out = graph_lower_smooth.SmoothLowess(graph_lower, "", 1.0)

    for i, index in enumerate(np.flatnonzero(nominal_y)):
        final_upper = np.sign(graph_upper_out.GetPointY(i)-1)*limit+1 if np.abs(graph_upper_out.GetPointY(i)-1) > limit else graph_upper_out.GetPointY(i)
        final_lower = np.sign(graph_lower_out.GetPointY(i)-1)*limit+1 if np.abs(graph_lower_out.GetPointY(i)-1) > limit else graph_lower_out.GetPointY(i)
        upper.SetBinContent(int(index)+1, final_upper*nominal_y[index])
        lower.SetBinContent(int(index)+1, final_lower*nominal_y[index])

with open(args.yamlpath, "r") as yamlfile:
    yaml_spec = yaml.safe_load(yamlfile)

year_to_plot = yaml_spec["year"]
lumi_unit_to_plot = yaml_spec["lumiunit"]
lumi_to_plot = yaml_spec["lumi"]
# Luminosity must be fb. If in pb, convert to fb first.
#if lumi_unit_to_plot == "pb": lumi_to_plot = lumi_to_plot*1000

analysis_name = "."
if "analysisname" in yaml_spec.keys(): analysis_name = yaml_spec["analysisname"]
if analysis_name != ".": os.system(f"mkdir {analysis_name}")

categories_to_plot = yaml_spec["categories"].keys()
treename_to_plot = yaml_spec["treename"]
basecut_to_plot = yaml_spec["basecut"]

genweight_to_plot = yaml_spec["genweight"]

mass_variables = [(v["variable"], v["bins"], v["range"]) for v in yaml_spec["distribution"]["variables"]]
for mass_variable, _, _ in mass_variables: os.system(f"mkdir {analysis_name}/{mass_variable}")

event_categories = [(e["name"], e["rule"]) for e in yaml_spec["distribution"]["event_categories"]]

tagger_name = yaml_spec["tagger"]["name"]
if "cutrule" in yaml_spec["tagger"].keys():
    tagger_cut_pass = yaml_spec["tagger"]["cutrule"]
    tagger_cut_fail = f"!({tagger_cut_pass})"
else:
    tagger_varname = yaml_spec["tagger"]["varname"]
    tagger_cut_pass = f'{tagger_varname}>={yaml_spec["tagger"]["cut"]}'
    tagger_cut_fail = f'{tagger_varname}<{yaml_spec["tagger"]["cut"]}'

unc_to_plot = [unc for unc in yaml_spec["uncertainties"].keys() if yaml_spec["uncertainties"][unc]["mode"] in ["factor", "file"]]

def extract_hist_dict(filelist, process, weight, mass_variable, mass_bins, mass_range, uncname="nominal"):
    cache_dict = {}
    
    for filecount, filepath in enumerate(filelist):
        print(f"Debug: filepath = {filepath}")
        cache_dict[filepath] = {}
        for event_catname, event_catrule in event_categories:
            if args.debug:
                print(f"Debug: event cat. name = {event_catname}")
                print(f"Debug: event cat. rule = {event_catrule}")
            cache_dict[filepath][event_catname] = {}
            for cat in categories_to_plot:
                if args.debug: print(f"Debug: cat = {cat}")
                if process not in yaml_spec["categories"][cat]["processes"]: continue
                category_cut = yaml_spec["categories"][cat]["cut"]
                
                histname = f"{process}_{mass_variable}_{uncname}_{filecount}_{event_catname}_{cat}"
                
                cache_dict[filepath][event_catname][cat] = {}
                cache_dict[filepath][event_catname][cat]["pass"] = extract_histogram(
                    filename=filepath, treename=treename_to_plot, var=mass_variable,
                    cut=f"({basecut_to_plot})&&({category_cut})&&({event_catrule})&&({tagger_cut_pass})",
                    weight=weight,
                    histname=histname+"_pass",
                    xbins=mass_bins, xmin=mass_range[0], xmax=mass_range[1]
                )
                cache_dict[filepath][event_catname][cat]["fail"] = extract_histogram(
                    filename=filepath, treename=treename_to_plot, var=mass_variable,
                    cut=f"({basecut_to_plot})&&({category_cut})&&({event_catrule})&&({tagger_cut_fail})",
                    weight=weight,
                    histname=histname+"_fail",
                    xbins=mass_bins, xmin=mass_range[0], xmax=mass_range[1]
                )
    return cache_dict

if "perfileweights" in yaml_spec.keys():
    for weightset in yaml_spec["perfileweights"]:
        branchname = weightset["name"]
        branchvalue = float(weightset["value"])
        for filename in weightset["files"]:
            rootfile = pyr.TFile(filename, "UPDATE")
            roottree = rootfile.Get(treename_to_plot)
            array_value = array('f', [0])
            new_branch = roottree.Branch(branchname, array_value, branchname + "/F")
            array_value[0] = branchvalue
            for _ in range(roottree.GetEntries()): new_branch.Fill()
            new_branch.ResetAddress()
            roottree.Write(treename_to_plot, pyr.TObject.kOverwrite)
            rootfile.Close()

hist_plots_per_processes_and_files = {}
for mass_variable, mass_bins, mass_range in mass_variables:
    hist_plots_per_processes_and_files[mass_variable] = {}
    for process in yaml_spec["processes"].keys():
        hist_plots_per_processes_and_files[mass_variable][process] = {}
        if process == "data": 
            print("Debug: Data")
            for filecount, filepath in enumerate(yaml_spec["processes"]["data"]["nominal_files"]):
                hist_plots_per_processes_and_files[mass_variable]["data"][filepath] = {}
                for event_catname, event_catrule in event_categories:
                    hist_plots_per_processes_and_files[mass_variable]["data"][filepath][event_catname] = {}
                    hist_plots_per_processes_and_files[mass_variable]["data"][filepath][event_catname]["pass"] = extract_histogram(
                        filename=filepath, treename=treename_to_plot, var=mass_variable,
                        cut=f"({basecut_to_plot})&&({event_catrule})&&({tagger_cut_pass})",
                        weight="1.",
                        histname=f"data_{mass_variable}_{filecount}_{event_catname}_pass",
                        xbins=mass_bins, xmin=mass_range[0], xmax=mass_range[1]
                    )
                    hist_plots_per_processes_and_files[mass_variable]["data"][filepath][event_catname]["fail"] = extract_histogram(
                        filename=filepath, treename=treename_to_plot, var=mass_variable,
                        cut=f"({basecut_to_plot})&&({event_catrule})&&({tagger_cut_fail})",
                        weight="1.",
                        histname=f"data_{mass_variable}_{filecount}_{event_catname}_fail",
                        xbins=mass_bins, xmin=mass_range[0], xmax=mass_range[1]
                    )
            continue
        if args.debug: print(f"Debug: process = {process}")
        
        hist_plots_per_processes_and_files[mass_variable][process]["nominal"] = {}
        for unc in unc_to_plot:
            if args.debug: print(f"Debug: unc = {unc}")
            hist_plots_per_processes_and_files[mass_variable][process][unc+"_up"] = {}
            hist_plots_per_processes_and_files[mass_variable][process][unc+"_down"] = {}
        
        weight_nominal = str(lumi_to_plot) + "*" + genweight_to_plot
        if "additional_weights" in yaml_spec["processes"][process].keys(): 
            weight_nominal += "*" + yaml_spec["processes"][process]["additional_weights"]
        if args.debug: print(weight_nominal)
        hist_plots_per_processes_and_files[mass_variable][process]["nominal"] = extract_hist_dict(
            yaml_spec["processes"][process]["nominal_files"],
            process=process, weight=weight_nominal, 
            mass_variable=mass_variable, mass_bins=mass_bins, mass_range=mass_range
        )
        if args.debug: print(hist_plots_per_processes_and_files[mass_variable][process]["nominal"])
        
        for unc in unc_to_plot:
            if args.debug: print(f"Debug: unc = {unc}")
            if yaml_spec["uncertainties"][unc]["mode"] == "factor":
                weight_uncup   = weight_nominal + "*" + yaml_spec["uncertainties"][unc]["up"]
                weight_uncdown = weight_nominal + "*" + yaml_spec["uncertainties"][unc]["down"]
                if args.debug: 
                    print(weight_uncup)
                    print(weight_uncdown)
                hist_plots_per_processes_and_files[mass_variable][process][unc+"_up"] = extract_hist_dict(
                    yaml_spec["processes"][process]["nominal_files"],
                    process=process, weight=weight_uncup, uncname=unc+"_up",
                    mass_variable=mass_variable, mass_bins=mass_bins, mass_range=mass_range
                )
                hist_plots_per_processes_and_files[mass_variable][process][unc+"_down"] = extract_hist_dict(
                    yaml_spec["processes"][process]["nominal_files"],
                    process=process, weight=weight_uncdown, uncname=unc+"_down",
                    mass_variable=mass_variable, mass_bins=mass_bins, mass_range=mass_range
                )
            elif yaml_spec["uncertainties"][unc]["mode"] == "file":
                hist_plots_per_processes_and_files[mass_variable][process][unc+"_up"] = extract_hist_dict(
                    yaml_spec["processes"][process]["unc_files"][unc]["up"],
                    process=process, weight=weight_nominal, uncname=unc+"_up",
                    mass_variable=mass_variable, mass_bins=mass_bins, mass_range=mass_range
                )
                hist_plots_per_processes_and_files[mass_variable][process][unc+"_down"] = extract_hist_dict(
                    yaml_spec["processes"][process]["unc_files"][unc]["down"],
                    process=process, weight=weight_nominal, uncname=unc+"_down",
                    mass_variable=mass_variable, mass_bins=mass_bins, mass_range=mass_range
                )

if args.diagnosis:
    for event_catname, event_catrule in event_categories:
        diagnosis_file = pyr.TFile(f"{analysis_name}/diagnosis_{event_catname}.root", "RECREATE")
        for mass_variable in hist_plots_per_processes_and_files.keys():
            for process in hist_plots_per_processes_and_files[mass_variable].keys():
                for uncvariant in hist_plots_per_processes_and_files[mass_variable][process].keys():
                    for filepath in hist_plots_per_processes_and_files[mass_variable][process][uncvariant].keys():
                        if event_catname in hist_plots_per_processes_and_files[mass_variable][process][uncvariant][filepath].keys():
                            for category in hist_plots_per_processes_and_files[mass_variable][process][uncvariant][filepath][event_catname].keys():
                                hist_plots_per_processes_and_files[mass_variable][process][uncvariant][filepath][event_catname][category]["pass"].Write()
                                hist_plots_per_processes_and_files[mass_variable][process][uncvariant][filepath][event_catname][category]["fail"].Write()
        diagnosis_file.Close()

hist_data_per_ptrange = {}
for mass_variable, mass_bins, mass_range in mass_variables:
    hist_data_per_ptrange[mass_variable] = {}
    for event_catname, event_catrule in event_categories:
        hist_data_per_ptrange[mass_variable][event_catname] = {"pass": {}, "fail": {}}
        hist_data_per_ptrange[mass_variable][event_catname]["pass"] = combine_histograms(
            [hist_plots_per_processes_and_files[mass_variable]["data"][filepath][event_catname]["pass"] for filepath in hist_plots_per_processes_and_files[mass_variable]["data"].keys()],
            finalname=f"data_{mass_variable}_{event_catname}_pass", 
            xbins=mass_bins, xmin=mass_range[0], xmax=mass_range[1]
        )
        hist_data_per_ptrange[mass_variable][event_catname]["fail"] = combine_histograms(
            [hist_plots_per_processes_and_files[mass_variable]["data"][filepath][event_catname]["fail"] for filepath in hist_plots_per_processes_and_files[mass_variable]["data"].keys()],
            finalname=f"data_{mass_variable}_{event_catname}_fail", 
            xbins=mass_bins, xmin=mass_range[0], xmax=mass_range[1]
        )

for mass_variable, mass_bins, mass_range in mass_variables:
    hist_plots_per_category = {}
    if args.debug: print(f"Debug: mass_variable = {mass_variable}")
    for category in categories_to_plot:
        if args.debug: print(f"Debug: category = {category}")
        hist_plots_per_category[category] = {}
        hist_plots_per_category[category]["nominal"] = {}
        for unc in unc_to_plot:
            hist_plots_per_category[category][unc+"_up"] = {}
            hist_plots_per_category[category][unc+"_down"] = {}
        
        for event_catname, event_catrule in event_categories:
            hist_plots_per_category[category]["nominal"][event_catname] = {}
            passing_list = []
            failing_list = []
            for process in yaml_spec["categories"][category]["processes"]:
                for filepath in hist_plots_per_processes_and_files[mass_variable][process]["nominal"].keys():
                    passing_list.append(hist_plots_per_processes_and_files[mass_variable][process]["nominal"][filepath][event_catname][category]["pass"])
                    failing_list.append(hist_plots_per_processes_and_files[mass_variable][process]["nominal"][filepath][event_catname][category]["fail"])
            hist_plots_per_category[category]["nominal"][event_catname]["pass"] = combine_histograms(
                histlist=passing_list,
                finalname=f"{category}_{mass_variable}_{event_catname}_pass_nominal",
                xbins=mass_bins, xmin=mass_range[0], xmax=mass_range[1]
            )
            hist_plots_per_category[category]["nominal"][event_catname]["fail"] = combine_histograms(
                histlist=failing_list,
                finalname=f"{category}_{mass_variable}_{event_catname}_fail_nominal",
                xbins=mass_bins, xmin=mass_range[0], xmax=mass_range[1]
            )
            
            for unc in unc_to_plot:
                hist_plots_per_category[category][unc+"_up"][event_catname] = {}
                hist_plots_per_category[category][unc+"_down"][event_catname] = {}
                passing_list_up   = []
                passing_list_down = []
                failing_list_up   = []
                failing_list_down = []
                for process in yaml_spec["categories"][category]["processes"]:
                    for filepath in hist_plots_per_processes_and_files[mass_variable][process][unc+"_up"].keys():
                        passing_list_up.append(hist_plots_per_processes_and_files[mass_variable][process][unc+"_up"][filepath][event_catname][category]["pass"])
                        failing_list_up.append(hist_plots_per_processes_and_files[mass_variable][process][unc+"_up"][filepath][event_catname][category]["fail"])
                    for filepath in hist_plots_per_processes_and_files[mass_variable][process][unc+"_down"].keys():
                        passing_list_down.append(hist_plots_per_processes_and_files[mass_variable][process][unc+"_down"][filepath][event_catname][category]["pass"])
                        failing_list_down.append(hist_plots_per_processes_and_files[mass_variable][process][unc+"_down"][filepath][event_catname][category]["fail"])
                        
                hist_plots_per_category[category][unc+"_up"][event_catname]["pass"] = combine_histograms(
                    histlist=passing_list_up,
                    finalname=f"{category}_{mass_variable}_{event_catname}_pass_{unc}Up",
                    xbins=mass_bins, xmin=mass_range[0], xmax=mass_range[1]
                )
                hist_plots_per_category[category][unc+"_up"][event_catname]["fail"] = combine_histograms(
                    histlist=failing_list_up,
                    finalname=f"{category}_{mass_variable}_{event_catname}_fail_{unc}Up",
                    xbins=mass_bins, xmin=mass_range[0], xmax=mass_range[1]
                )
                
                hist_plots_per_category[category][unc+"_down"][event_catname]["pass"] = combine_histograms(
                    histlist=passing_list_down,
                    finalname=f"{category}_{mass_variable}_{event_catname}_pass_{unc}Down",
                    xbins=mass_bins, xmin=mass_range[0], xmax=mass_range[1]
                )
                hist_plots_per_category[category][unc+"_down"][event_catname]["fail"] = combine_histograms(
                    histlist=failing_list_down,
                    finalname=f"{category}_{mass_variable}_{event_catname}_fail_{unc}Down",
                    xbins=mass_bins, xmin=mass_range[0], xmax=mass_range[1]
                )
    
    for category in hist_plots_per_category.keys():
        for unc in unc_to_plot:
            if "lowess" in yaml_spec["uncertainties"][unc].keys():
                for pt_range in hist_plots_per_category[category][unc+"_up"].keys():
                    for passorfail in hist_plots_per_category[category][unc+"_up"][pt_range].keys():
                        create_lowess(nominal=hist_plots_per_category[category]["nominal"][pt_range][passorfail], 
                                      upper=hist_plots_per_category[category][unc+"_up"][pt_range][passorfail], 
                                      lower=hist_plots_per_category[category][unc+"_down"][pt_range][passorfail], 
                                      symmetrize=yaml_spec["uncertainties"][unc]["lowess"]["symmetrize"],
                                      limit=yaml_spec["uncertainties"][unc]["lowess"]["limit"]
                                      )

    for category in hist_plots_per_category.keys():
        print("=====================")
        print(f"Category: {category}")
        for unc in hist_plots_per_category[category].keys():
            print(f"Uncertainty: {unc}")
            for pt_range in hist_plots_per_category[category][unc].keys():
                print(f"pT range: {pt_range}")
                for passorfail in hist_plots_per_category[category][unc][pt_range].keys():
                    print(passorfail)
                    print(hist_plots_per_category[category][unc][pt_range][passorfail])

    analysis_obj_collection = {}
    #for i, pt_range in enumerate(pt_ranges_to_plot):
    for event_catname, event_catrule in event_categories:
        analysis_hist_obj = AnalysisHistogram(categories_to_plot, mass_bins, mass_range[0], mass_range[1], treename_to_plot)
        #pt_range_name = pt_ranges_name[i]
        analysis_hist_obj.add_data_hist(hist=hist_data_per_ptrange[mass_variable][event_catname]["pass"], isPass=True)
        analysis_hist_obj.add_data_hist(hist=hist_data_per_ptrange[mass_variable][event_catname]["fail"], isPass=False)
        for category in categories_to_plot:
            analysis_hist_obj.add_nominal_hist(
                category=category, 
                histobj=hist_plots_per_category[category]["nominal"][event_catname]["pass"],
                isPass=True
            )
            analysis_hist_obj.add_nominal_hist(
                category=category, 
                histobj=hist_plots_per_category[category]["nominal"][event_catname]["fail"],
                isPass=False
            )
            for unc in unc_to_plot:
                analysis_hist_obj.define_unc(category=category, unc=unc)
                analysis_hist_obj.add_unc_hist(
                    category=category, unc=unc,
                    histobj=hist_plots_per_category[category][unc+"_up"][event_catname]["pass"], 
                    isUp=True, isPass=True
                )
                analysis_hist_obj.add_unc_hist(
                    category=category, unc=unc,
                    histobj=hist_plots_per_category[category][unc+"_down"][event_catname]["pass"], 
                    isUp=False, isPass=True
                )
                analysis_hist_obj.add_unc_hist(
                    category=category, unc=unc,
                    histobj=hist_plots_per_category[category][unc+"_up"][event_catname]["fail"], 
                    isUp=True, isPass=False
                )
                analysis_hist_obj.add_unc_hist(
                    category=category, unc=unc,
                    histobj=hist_plots_per_category[category][unc+"_down"][event_catname]["fail"], 
                    isUp=False, isPass=False
                )
        analysis_obj_collection[event_catname] = analysis_hist_obj
    if args.debug: print(analysis_obj_collection)
    for key in analysis_obj_collection.keys():
        if args.debug: print(analysis_obj_collection[key].__dict__)
        analysis_obj_collection[key].save_histograms(f"{analysis_name}/{mass_variable}/{key}.root")
        
        number_of_categories = len(analysis_obj_collection[key].categories)
        data_pass_count = analysis_obj_collection[key].data_hist["pass"].Integral(1, analysis_obj_collection[key].data_hist["pass"].GetNbinsX())
        data_fail_count = analysis_obj_collection[key].data_hist["fail"].Integral(1, analysis_obj_collection[key].data_hist["fail"].GetNbinsX())
        mc_pass_count = 0
        mc_fail_count = 0
        for category in analysis_obj_collection[key].categories:
            pass_bins = analysis_obj_collection[key].nom_hist[category]["pass"].GetNbinsX()
            fail_bins = analysis_obj_collection[key].nom_hist[category]["fail"].GetNbinsX()
            mc_pass_count += analysis_obj_collection[key].nom_hist[category]["pass"].Integral(1, pass_bins)
            mc_fail_count += analysis_obj_collection[key].nom_hist[category]["fail"].Integral(1, fail_bins)
        norm_match_mc_to_data = (data_pass_count + data_fail_count) / (mc_pass_count + mc_fail_count)
        
        print(f"Printing datacard {key}.txt")
        with open(f"{analysis_name}/{mass_variable}/{key}.txt", 'w') as combine_card_file:
            combine_lines = []
            
            combine_lines.append("imax 2 (two channels, pass and fail)\n")
            combine_lines.append(f"jmax {number_of_categories-1} ({number_of_categories} categories minus 1)\n")
            combine_lines.append("kmax * (automatic number of nuisance parameters)\n")
            combine_lines.append("----------\n")
            
            combine_lines.append(f"shapes data_obs pass {key}.root data_{mass_variable}_{key}_pass\n")
            combine_lines.append(f"shapes * pass {key}.root $PROCESS_{mass_variable}_{key}_pass_nominal $PROCESS_{mass_variable}_{key}_pass_$SYSTEMATIC\n")
            combine_lines.append(f"shapes data_obs fail {key}.root data_{mass_variable}_{key}_fail\n")
            combine_lines.append(f"shapes * fail {key}.root $PROCESS_{mass_variable}_{key}_fail_nominal $PROCESS_{mass_variable}_{key}_fail_$SYSTEMATIC\n")
            combine_lines.append("----------\n")
            
            combine_lines.append("bin\tpass\tfail\n")
            combine_lines.append(f"observation\t{data_pass_count:.0f}\t{data_fail_count:.0f}\n")
            combine_lines.append("----------\n")
            
            combine_lines.append("# automatic counting of MC events\n")
            combine_lines.append("bin\t" + "pass\t"*number_of_categories + "fail\t"*number_of_categories + '\n')
            combine_lines.append("process\t"+'\t'.join(list(analysis_obj_collection[key].categories)*2) + '\n')
            combine_lines.append("process\t"+'\t'.join(list(map(str, range(number_of_categories)))*2) + '\n')
            combine_lines.append("rate\t"+"-1\t"*number_of_categories*2 + '\n')
            combine_lines.append("----------\n")
            
            for unc in yaml_spec["uncertainties"].keys():
                unc_line = unc + '\t'
                unc_size = ""
                unc_size_number = yaml_spec["uncertainties"][unc]["size"] if "size" in yaml_spec["uncertainties"][unc] else 1
                if yaml_spec["uncertainties"][unc]["mode"] in ["factor", "file"]: 
                    unc_line += "shape\t"
                else: 
                    unc_line += yaml_spec["uncertainties"][unc]["mode"] + '\t'
                if "category" in yaml_spec["uncertainties"][unc].keys():
                    for category in analysis_obj_collection[key].categories:
                        if category == yaml_spec["uncertainties"][unc]["category"]: unc_size += f"{unc_size_number}\t"
                        else: unc_size += "-\t"
                    unc_size = unc_size + unc_size
                else:
                    unc_size = '\t'.join([f"{unc_size_number}"]*number_of_categories*2)
                combine_lines.append(unc_line+unc_size+'\n')
            combine_lines.append("# normalisation factor to match MC and data\n")
            combine_lines.append("# freezes automatically\n")
            combine_lines.append(f"norm_match_mc_data rateParam * * {norm_match_mc_to_data:.6f}\n")
            combine_lines.append("nuisance edit freeze norm_match_mc_data\n")
            combine_lines.append("\n")
            combine_lines.append("# activating autoMCStats\n")
            combine_lines.append("* autoMCStats 0\n")
            
            combine_card_file.writelines(combine_lines)

    print(f"Printing combine script file combine_script.sh")
    with open(f"{analysis_name}/{mass_variable}/combine_script.sh", 'w') as combine_script_file:
        combine_script_file.write("#!/bin/bash\n")
        combine_script_file.write("# Converting datacards to workspace file for portability :-)\n")
        for key in analysis_obj_collection.keys():
            combine_script_file.write(f"text2workspace.py -m 125 -P HiggsAnalysis.CombinedLimit.TagAndProbeExtended:tagAndProbe {key}.txt -o workspace_{key}.root --PO=categories={','.join(categories_to_plot)}\n")

print("===========================")
print("All done! :-)")
print("See COMBINE_README.md for more info on combine script usage.")