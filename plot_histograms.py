import argparse
import numpy as np
import ROOT as pyr
from matplotlib import pyplot as plt
import yaml
import mplhep as hep

plt.style.use(hep.style.CMS)

def hist_to_bins(hist):
    return [hist.GetBinLowEdge(i) for i in range(1, hist.GetNbinsX()+2)]

def hist_to_array(hist, isData=False):
    res = []
    err = []
    for i in range(1, hist.GetNbinsX()+1):
        res.append(hist.GetBinContent(i))
        if isData: err.append((hist.GetBinErrorUp(i), hist.GetBinErrorLow(i)))
        else: err.append(hist.GetBinError(i))
    res = np.array(res)
    err = np.array(err).transpose()
    return res, err

def graph_to_array(graph):
    res = []
    err = []
    for i in range(graph.GetN()):
        res.append(graph.GetPointY(i))
        err.append((graph.GetErrorYhigh(i), graph.GetErrorYlow(i)))
    res = np.array(res)
    err = np.array(err).transpose()
    return res, err

parser = argparse.ArgumentParser()
parser.add_argument("yamlpath", help="input YAML file path, not the same as input for make_histograms.py")
args = parser.parse_args()

with open(args.yamlpath, "r") as yamlfile:
    yaml_spec = yaml.safe_load(yamlfile)

def plot_prefit(array_mc, array_mc_error, array_mc_sum, array_data, array_data_err, histbins, legendtitle, filename):
    fig = plt.figure(figsize=(12, 12), facecolor="white")
    main_ax = plt.subplot2grid((5, 1), (0, 0), rowspan=4)
    ratio_ax = plt.subplot2grid((5, 1), (4, 0))
    fig.subplots_adjust(hspace=0)
    
    baseline = np.zeros(array_mc_sum.shape)
    for category, category_config in yaml_spec["categories"].items():
        main_ax.stairs(
            array_mc[category] + baseline,
            histbins,
            baseline=baseline,
            fill=True,
            label=category_config["propername"],
            facecolor=category_config["color"]
        )
        main_ax.errorbar(
            x=[(histbins[i]+histbins[i+1])/2 for i in range(len(histbins)-1)],
            y=array_mc[category] + baseline,
            yerr=array_mc_error[category],
            fmt="none",
            ecolor=category_config["color"],
            capsize=5
        )
        baseline += array_mc[category]
    main_ax.errorbar(
        x=[(histbins[i]+histbins[i+1])/2 for i in range(len(histbins)-1)],
        y=array_data,
        yerr=array_data_err,
        fmt="o",
        color="black",
        label="Data",
        capsize=5
    )
    
    handles, labels = main_ax.get_legend_handles_labels()
    order = [len(labels)-1] + list(range(len(labels)-1))
    main_ax.legend(
        [handles[idx] for idx in order],
        [labels[idx] for idx in order],
        title=legendtitle
    )
    
    ratio_ax.hlines(1, min(histbins), max(histbins), linestyles="--", colors="#9c9ca1")
    ratio_ax.errorbar(
        x=[(histbins[i]+histbins[i+1])/2 for i in range(len(histbins)-1)],
        y=array_data/array_mc_sum,
        yerr=array_data_err/array_mc_sum,
        fmt="o",
        color="black",
        capsize=5
    )
    
    main_ax.set_xlim((min(histbins), max(histbins)))
    main_ax.set_ylim((0, max(array_mc_sum)*1.5))
    main_ax.set_xticklabels([])
    main_ax.yaxis.get_major_ticks()[0].label1.set_visible(False)
    main_ax.set_ylabel("Events")
    
    ratio_ax.set_xlim((min(histbins), max(histbins)))
    ratio_ax.set_ylim((0.75, 1.25))
    ratio_ax.set_yticks((0.75, 1, 1.25))
    ratio_ax.set_xlabel(yaml_spec["xlabel"])
    ratio_ax.set_ylabel("Data/MC")
    
    hep.cms.label(
        llabel="Preliminary",
        lumi=yaml_spec["lumi"],
        ax=main_ax,
    )
    
    fig.savefig(filename, bbox_inches="tight")

def plot_postfit(array_prefit_mc, array_prefit_mc_error, array_prefit_mc_sum, array_prefit_mc_sum_error, array_postfit_mc, array_postfit_mc_error, array_postfit_mc_sum, array_postfit_mc_sum_error, array_data, array_data_err, histbins, legendtitle, filename):
    fig = plt.figure(figsize=(12, 12), facecolor="white")
    main_ax = plt.subplot2grid((5, 1), (0, 0), rowspan=4)
    ratio_ax = plt.subplot2grid((5, 1), (4, 0))
    fig.subplots_adjust(hspace=0)
    
    for category, category_config in yaml_spec["categories"].items():
        main_ax.stairs(
            array_prefit_mc[category],
            histbins,
            fill=False,
            color=category_config["color"],
            linestyle=":",
            linewidth=2
        )
        prefit_errbar = main_ax.errorbar(
            x=[(histbins[i]+histbins[i+1])/2 for i in range(len(histbins)-1)],
            y=array_prefit_mc[category],
            yerr=array_prefit_mc_error[category],
            fmt="none",
            ecolor=category_config["color"],
            linestyle=":",
            linewidth=2
        )
        prefit_errbar[-1][0].set_linestyle(":")
        main_ax.stairs(
            array_postfit_mc[category],
            histbins,
            fill=False,
            label=category_config["propername"],
            color=category_config["color"],
            linestyle="-",
            linewidth=2
        )
        main_ax.errorbar(
            x=[(histbins[i]+histbins[i+1])/2 for i in range(len(histbins)-1)],
            y=array_postfit_mc[category],
            yerr=array_postfit_mc_error[category],
            fmt="none",
            ecolor=category_config["color"],
            linestyle="-",
            capsize=5,
            linewidth=2
        )
    
    main_ax.stairs(
        array_prefit_mc_sum,
        histbins,
        fill=False,
        color="#717581",
        label="Prefit",
        linestyle=":",
        linewidth=2
    )
    total_prefit_errbar = main_ax.errorbar(
        x=[(histbins[i]+histbins[i+1])/2 for i in range(len(histbins)-1)],
        y=array_prefit_mc_sum,
        yerr=array_prefit_mc_sum_error,
        fmt="none",
        ecolor="#717581",
        linestyle=":",
        linewidth=2
    )
    total_prefit_errbar[-1][0].set_linestyle(":")
    main_ax.stairs(
        array_postfit_mc_sum,
        histbins,
        fill=False,
        label="Total SM",
        color="#717581",
        linestyle="-",
        linewidth=2
    )
    main_ax.stairs(
        array_postfit_mc_sum+array_postfit_mc_sum_error[0],
        histbins,
        baseline=array_postfit_mc_sum-array_postfit_mc_sum_error[1],
        fill=False,
        label="Postfit unc.",
        color="#717581",
        hatch="///",
        linewidth=0
    )
    
    main_ax.errorbar(
        x=[(histbins[i]+histbins[i+1])/2 for i in range(len(histbins)-1)],
        y=array_data,
        yerr=array_data_err,
        fmt="o",
        color="black",
        label="Data",
        capsize=5,
        linewidth=2
    )
    
    handles, labels = main_ax.get_legend_handles_labels()
    order = [len(labels)-1] + list(range(len(labels)-4)) + [len(labels)-3, len(labels)-2, len(labels)-4]
    main_ax.legend(
        [handles[idx] for idx in order],
        [labels[idx] for idx in order],
        title=legendtitle
    )
    
    ratio_ax.hlines(1, min(histbins), max(histbins), linestyles="--", colors="#9c9ca1")
    ratio_ax.stairs(
        (array_postfit_mc_sum+array_postfit_mc_sum_error[0])/array_postfit_mc_sum,
        histbins,
        baseline=(array_postfit_mc_sum-array_postfit_mc_sum_error[1])/array_postfit_mc_sum,
        fill=False,
        label="Postfit unc.",
        color="black",
        hatch="///",
        linewidth=0
    )
    ratio_ax.errorbar(
        x=[(histbins[i]+histbins[i+1])/2 for i in range(len(histbins)-1)],
        y=array_data/array_postfit_mc_sum,
        yerr=array_data_err/array_postfit_mc_sum,
        fmt="o",
        color="black",
        capsize=5
    )
    
    main_ax.set_xlim((min(histbins), max(histbins)))
    main_ax.set_ylim((0, max(array_postfit_mc_sum)*1.5))
    main_ax.set_xticklabels([])
    main_ax.yaxis.get_major_ticks()[0].label1.set_visible(False)
    main_ax.set_ylabel("Events")
    
    ratio_ax.set_xlim((min(histbins), max(histbins)))
    ratio_ax.set_ylim((0.75, 1.25))
    ratio_ax.set_yticks((0.75, 1, 1.25))
    ratio_ax.set_xlabel(yaml_spec["xlabel"])
    ratio_ax.set_ylabel("Data/MC")
    
    hep.cms.label(
        llabel="Preliminary",
        lumi=yaml_spec["lumi"],
        ax=main_ax,
    )
    
    fig.savefig(filename, bbox_inches="tight")

for ptrange_dict in yaml_spec["ptranges"]:
    ptrange_name = ptrange_dict["name"]
    if "propername" in ptrange_dict.keys(): ptrange_propername = ptrange_dict["propername"]
    else: ptrange_propername = ptrange_name
    prefitfile  = pyr.TFile(ptrange_dict["prefitfile"])
    postfitfile = pyr.TFile(ptrange_dict["postfitfile"])
    
    hist_prefit_data_pass = prefitfile.Get(f"data_{ptrange_name}_pass")
    hist_prefit_data_fail = prefitfile.Get(f"data_{ptrange_name}_fail")
    hist_prefit_mc_pass = {}
    hist_prefit_mc_fail = {}
    for category in yaml_spec["categories"].keys():
        hist_prefit_mc_pass[category] = prefitfile.Get(f"{category}_{ptrange_name}_pass_nominal")
        hist_prefit_mc_fail[category] = prefitfile.Get(f"{category}_{ptrange_name}_fail_nominal")
    hist_prefit_mc_pass_sum = hist_prefit_mc_pass[list(yaml_spec["categories"].keys())[0]].Clone(f"total_{ptrange_name}_pass_nominal")
    hist_prefit_mc_fail_sum = hist_prefit_mc_fail[list(yaml_spec["categories"].keys())[0]].Clone(f"total_{ptrange_name}_fail_nominal")
    hist_prefit_mc_pass_sum.Reset("ICES")
    hist_prefit_mc_fail_sum.Reset("ICES")
    for category in yaml_spec["categories"].keys():
        hist_prefit_mc_pass_sum.Add(hist_prefit_mc_pass[category])
        hist_prefit_mc_fail_sum.Add(hist_prefit_mc_fail[category])
    
    hist_postfit_data_pass = postfitfile.Get(f"shapes_prefit/pass/data")
    hist_postfit_data_fail = postfitfile.Get(f"shapes_prefit/fail/data")
    hist_postfit_mc_pass = {}
    hist_postfit_mc_fail = {}
    hist_postfit_mc_pass_prefit = {}
    hist_postfit_mc_fail_prefit = {}
    for category in yaml_spec["categories"].keys():
        hist_postfit_mc_pass[category] = postfitfile.Get(f"shapes_fit_s/pass/{category}")
        hist_postfit_mc_fail[category] = postfitfile.Get(f"shapes_fit_s/fail/{category}")
        hist_postfit_mc_pass_prefit[category] = postfitfile.Get(f"shapes_prefit/pass/{category}")
        hist_postfit_mc_fail_prefit[category] = postfitfile.Get(f"shapes_prefit/fail/{category}")
    hist_postfit_mc_pass_sum = postfitfile.Get("shapes_fit_s/pass/total")
    hist_postfit_mc_fail_sum = postfitfile.Get("shapes_fit_s/fail/total")
    hist_postfit_mc_pass_prefit_sum = postfitfile.Get("shapes_prefit/pass/total")
    hist_postfit_mc_fail_prefit_sum = postfitfile.Get("shapes_prefit/fail/total")
    
    array_prefit_data_pass, array_prefit_data_pass_err = hist_to_array(hist_prefit_data_pass, isData=True)
    array_prefit_data_fail, array_prefit_data_fail_err = hist_to_array(hist_prefit_data_fail, isData=True)
    array_prefit_mc_pass = {}
    array_prefit_mc_fail = {}
    array_prefit_mc_pass_error = {}
    array_prefit_mc_fail_error = {}
    for category in yaml_spec["categories"].keys():
        array_prefit_mc_pass[category], array_prefit_mc_pass_error[category] = hist_to_array(hist_prefit_mc_pass[category])
        array_prefit_mc_fail[category], array_prefit_mc_fail_error[category] = hist_to_array(hist_prefit_mc_fail[category])
    array_prefit_mc_pass_sum, array_prefit_mc_pass_sum_err = hist_to_array(hist_prefit_mc_pass_sum)
    array_prefit_mc_fail_sum, array_prefit_mc_fail_sum_err = hist_to_array(hist_prefit_mc_fail_sum)
    
    array_postfit_data_pass, array_postfit_data_pass_err = graph_to_array(hist_postfit_data_pass)
    array_postfit_data_fail, array_postfit_data_fail_err = graph_to_array(hist_postfit_data_fail)
    array_postfit_mc_pass = {}
    array_postfit_mc_fail = {}
    array_postfit_mc_pass_err = {}
    array_postfit_mc_fail_err = {}
    array_postfit_mc_pass_prefit = {}
    array_postfit_mc_fail_prefit = {}
    array_postfit_mc_pass_prefit_err = {}
    array_postfit_mc_fail_prefit_err = {}
    for category in yaml_spec["categories"].keys():
        array_postfit_mc_pass[category], array_postfit_mc_pass_err[category] = hist_to_array(hist_postfit_mc_pass[category])
        array_postfit_mc_fail[category], array_postfit_mc_fail_err[category] = hist_to_array(hist_postfit_mc_fail[category])
        array_postfit_mc_pass_prefit[category], array_postfit_mc_pass_prefit_err[category] = hist_to_array(hist_postfit_mc_pass_prefit[category])
        array_postfit_mc_fail_prefit[category], array_postfit_mc_fail_prefit_err[category] = hist_to_array(hist_postfit_mc_fail_prefit[category])
    array_postfit_mc_pass_sum, array_postfit_mc_pass_sum_err = hist_to_array(hist_postfit_mc_pass_sum)
    array_postfit_mc_fail_sum, array_postfit_mc_fail_sum_err = hist_to_array(hist_postfit_mc_fail_sum)
    array_postfit_mc_pass_prefit_sum, array_postfit_mc_pass_prefit_sum_err = hist_to_array(hist_postfit_mc_pass_prefit_sum)
    array_postfit_mc_fail_prefit_sum, array_postfit_mc_fail_prefit_sum_err = hist_to_array(hist_postfit_mc_fail_prefit_sum)
    
    histbins_prefit_pass = hist_to_bins(hist_prefit_mc_pass_sum)
    histbins_prefit_fail = hist_to_bins(hist_prefit_mc_fail_sum)
    histbins_postfit_pass = hist_to_bins(hist_postfit_mc_pass_sum)
    histbins_postfit_fail = hist_to_bins(hist_postfit_mc_fail_sum)
    
    #print("---")
    #print(array_prefit_data_pass)
    #print("err", array_prefit_data_pass_err)
    #print(array_prefit_data_fail)
    #print("err", array_prefit_data_fail_err)
    #for category in yaml_spec["categories"].keys():
    #    print(category)
    #    print(array_postfit_mc_pass[category])
    #    print("err", array_postfit_mc_pass_err[category])
    #    print(array_postfit_mc_fail[category])
    #    print("err", array_postfit_mc_fail_err[category])
    #    print(array_postfit_mc_pass_prefit[category])
    #    print("err", array_postfit_mc_pass_prefit_err[category])
    #    print(array_postfit_mc_fail_prefit[category])
    #    print("err", array_postfit_mc_fail_prefit_err[category])
    #print(array_postfit_mc_pass_sum)
    #print("err", array_postfit_mc_pass_sum_err)
    #print(array_postfit_mc_fail_sum)
    #print("err", array_postfit_mc_fail_sum_err)
    #print(array_postfit_mc_pass_prefit_sum)
    #print("err", array_postfit_mc_pass_prefit_sum_err)
    #print(array_postfit_mc_fail_prefit_sum)
    #print("err", array_postfit_mc_fail_prefit_sum_err)
    
    # First plot prefit
    plot_prefit(
        array_prefit_mc_pass, 
        array_prefit_mc_pass_error, 
        array_prefit_mc_pass_sum, 
        array_prefit_data_pass, 
        array_prefit_data_pass_err, 
        histbins_prefit_pass, 
        ptrange_propername + ", pass", 
        f"{yaml_spec['savedir']}/prefit_pass_{ptrange_name}.png"
    )
    plot_prefit(
        array_prefit_mc_fail, 
        array_prefit_mc_fail_error, 
        array_prefit_mc_fail_sum, 
        array_prefit_data_fail, 
        array_prefit_data_fail_err, 
        histbins_prefit_fail, 
        ptrange_propername + ", fail", 
        f"{yaml_spec['savedir']}/prefit_fail_{ptrange_name}.png"
    )
    plot_postfit(
        array_postfit_mc_pass_prefit, 
        array_postfit_mc_pass_prefit_err, 
        array_postfit_mc_pass_prefit_sum, 
        array_postfit_mc_pass_prefit_sum_err, 
        array_postfit_mc_pass, 
        array_postfit_mc_pass_err, 
        array_postfit_mc_pass_sum, 
        array_postfit_mc_pass_sum_err, 
        array_postfit_data_pass, 
        array_postfit_data_pass_err, 
        histbins_prefit_pass,
        ptrange_propername + ", pass", 
        f"{yaml_spec['savedir']}/postfit_pass_{ptrange_name}.png"
    )
    plot_postfit(
        array_postfit_mc_fail_prefit, 
        array_postfit_mc_fail_prefit_err, 
        array_postfit_mc_fail_prefit_sum, 
        array_postfit_mc_fail_prefit_sum_err, 
        array_postfit_mc_fail, 
        array_postfit_mc_fail_err, 
        array_postfit_mc_fail_sum, 
        array_postfit_mc_fail_sum_err, 
        array_postfit_data_fail, 
        array_postfit_data_fail_err, 
        histbins_prefit_fail,
        ptrange_propername + ", fail", 
        f"{yaml_spec['savedir']}/postfit_fail_{ptrange_name}.png"
    )