# TopSF
A _better_ framework for scale factor measurement, designed with flexibility in mind.

**This repository is under construction!**

## Why use this?
- Supports _any_ ROOT ntuple structure. No need to conform to any hardcoded format.
- Supports _any_ tagger that gives one discriminator value, such as cut-based, BDT, or neural networks.
- Supports _any_ file name format, _any_ number of processes, and _any_ number of tagging categories you have.
- _Everything_ can be defined in one YAML file for datacard creation. No more surprises hidden deep in raw code.
- Also contains _helpful_ HiggsCombine script to give you an idea of what you can do with the output datacard.

## Requirements
- Python 3 with pyYAML, numpy, matplotlib, and mplhep
- ROOT with pyROOT interface

## Usage
This framework contains two main Python scripts:
```bash
# Make distributions for scale factor fitting.
python make_histogram.py YAML_FILE [--diagnosis]

# Make prefit and postfit plots.
python plot_histograms.py PLOT_YAML_FILE
```
`make_histogram.py` is the main script that generates the _final_ 1D distributions containing events passing and failing the designated tagger. It requires a YAML input file detailing everything regarding the setup, such as input ROOT file location, processes and tagging categories involved, and uncertainty definitions. 

The output files from this script are, per one pT category,
- ROOT file containing 1D distributions
- an accompanying combine card for that pT category
- if `--diagnosis` option is present, diagnosis ROOT file containing all distributions created from each input ROOT file _for MC_, arranged by the file name order

This script will also create one helpful bash script invoking `text2workspace` program, which can be used on machines with HiggsCombine set up.

Normally, to save time, other frameworks may generate the intermediate 2D histogram templates (containing jet pT versus jet mass distribution) for fast datacard generation in case the user wants to adjust the jet pT range. Unfortunately this may lead to bugs since the 2D histogram may not always have the exact pT ranges encoded. To avoid this surprise, **this script will only generate 1D distribution and no intermediate 2D histogram templates**.

The output files from this script should be used to measure scale factors using HiggsCombine. Refer to `COMBINE_README.md` file for more details.

Once you run the datacard, generated from the first Python file, using FitDiagnostics method in Higgs Combine, `plot_histograms.py` script can take the output file and generate both prefit and postfit plots at the same time.

### YAML input specifications
Refer to `specfile_test_2018.yaml` and `plotfile_2018.yaml` for examples of YAML file structures for `make_histograms.py` and `plot_histograms.py`.

#### `make_histograms.py`
The YAML input for `make_histograms.py` must contain details about your input ROOT files and the configuation you would like to have as follows:
- `year`: Analysis year. **Required.**
- `lumi`: Integrated luminosity. **Required.**
- `lumiunit`: Unit of integrated luminosity. **Required.**
- `genweight`: Generator weight definitions. Can be written as an expression. **Required.**
- `treename`: Tree name to look for in the input ROOT files. **Required, must have the same name for all ROOT files.**
- `processes`: List of processes to be added into tagging categories. **Required.** The structure must be as follows:
    - **One `data` key is required.** This specifies the input ROOT files for data. Another key `nominal_files` must be present in `data`, and must contain a list of input file paths.
    - Any number of processes, each as one key. Inside it must contain the keys `nominal_files` that specify the file paths for nominal event distributions, and `unc_files` that specify the file paths for event distributions for shape-based uncertainties, such as JES or JER.
- `basecut`: Base cut to be applied to all events. **Required.**
- `categories`: _Tagging categories_ to be added to the datacard. Each key represents a tagging category, and must contain the list of processes (under `processes` key) and further cuts for that category (under `cut` key).
- `tagger`: Details on the designated tagger
    - `name`: Name of the tagger to be used in files
    - `varname`: Name of the discrimnator _as seen in the ntuple files_.
    - `cut`: Discriminator cut for the tagger. Passing and failing events are defined as events with discriminator values higher and lower than this cut respectively.
- `distribution`: Details on the final distribution for scale factor measurement.
    - `mass_variable`: Target variable to be populated in the distribution. Usually jet mass is used.
    - `mass_range`: Range of the target variable.
    - `mass_bins`: Number of bins in the distribution.
    - `pt_variable`: Variable for pT value, used for different jet pT categories.
    - `pt_ranges`: Ranges of pt variables for datacards to be generated. Must be a list. If an element in the list is a list with two numbers, the jet pT category will have the default name as `"{pt[0]}to{pt[1]}"`. If an element in the list is a list with three elements, the first and second number will be used as pT range, while the third _string_ will be treated as jet pT category name. **The pT range name defined here will be used as output file names, such as `600to1200.root` and `600to1200.txt`.**
- `uncertainties`: Details on uncertainties to be added into the datacard. Each key is a uncertainty name, and should contain `mode` key inside. Currently `mode` supports `lnN`, `factor`, and `file`, and has different behaviours as follows:
    - `lnN`: Log-normal uncertainty. If no `category` key is present, this uncertainty will be applied to all _tagging categories_ with the specified `size` value. Specify tagging categories to apply this uncertainty using `category` key.
    - `factor`: Shape uncertainty calculated from _nominal_ input files with the designated expression. Must contain keys `up` and `down`.
    - `file`: Shape uncertainty calculated from files specified in `unc_files`. In this case, the uncertainty name must be the same as specified in `processes` key.

#### `plot_histograms.py`
The YAML input for `plot_histograms.py` has a different structure, aimed at plotting the prefit and postfit histograms as follows:
- `lumi`: Integrated luminosity. **Required.**
- `xlabel`: Label on the x-axis. **Required.** The range of x-axis is automatically determined by the histogram bins.
- `savedir`: Save directory **Required**
- `categories`: _Tagging categories_ as defined in the YAML input file for `make_histograms.py`. However, each key (name of the tagging category) must contain the following sub-keys:
    - `color`: Histogram colour, used in both prefit and postfit plots
    - `propername`: Name to be shown in the legend
- `ptranges`: _List_ of $p_T$ categories to be plotted. Must contain the following sub-keys:
    - `name`: $p_T$ category name as defined in `pt_ranges` option in the input YAML file for `make_histograms.py`
    - `propername`: Text to be added in the legend of the plot, _not required_
    - `prefitfile`: Path to prefit file, or the input ROOT file generated from `make_histograms.py`
    - `postfitfile`: Path to postfit file, or the output ROOT file from FitDiagnostics method of Higgs Combine

### What's inside the output ROOT file
The output ROOT file contains _all_ 1D histograms including passing and failing distributions. The naming convention is as follows:
- For data histograms, `data_{PT_CATEGORY}_pass` or `data_{PT_CATEGORY}_fail`
- For MC histograms, 
    - `{TAGGING_CATEGORY}_{PT_CATEGORY}_pass_nominal`
    - `{TAGGING_CATEGORY}_{PT_CATEGORY}_pass_{UNCERTAINTY}Up`
    - `{TAGGING_CATEGORY}_{PT_CATEGORY}_pass_{UNCERTAINTY}Down`
    - `{TAGGING_CATEGORY}_{PT_CATEGORY}_fail_nominal`
    - `{TAGGING_CATEGORY}_{PT_CATEGORY}_fail_{UNCERTAINTY}Up`
    - `{TAGGING_CATEGORY}_{PT_CATEGORY}_fail_{UNCERTAINTY}Down`

If `--diagnosis` option is turned on for `make_histograms.py`, another ROOT file, per pT category, will be created with the name `diagnosis_{PT_CATEGORY}.root` The naming convention in this file is 
`{PROCESS}_{UNCERTAINTY}_{up|down}_{FILE_INDEX}_{PT_CATEGORY}_{TAGGING_CATEGORY}_{pass|fail}`
where `FILE_INDEX` represents the order of the input file specified in the YAML file.