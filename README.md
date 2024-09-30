# TopSF
A _better_ framework for scale factor measurement, designed with flexibility in mind.

**For feature request or bugs please make an issue in this GitHub repository. Feedback is always appreciated here. Thanks!**

## Why use this?
- Supports _any_ ROOT ntuple structure. No need to conform to any hardcoded format.
- Supports _any_ tagger that gives one discriminator value, such as cut-based, BDT, or neural networks.
- Supports _any_ file name format, _any_ number of processes, and _any_ number of tagging categories you have.
- Supports _any_ kind of event category definition. You are not limited to pT ranges. You can define event categories in _any_ way with _any_ number of variables.
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

The output files from this script are, per one event category,
- ROOT file containing 1D distributions
- an accompanying combine card for that event category
- if `--diagnosis` option is present, diagnosis ROOT file containing all distributions created from each input ROOT file _for MC_, arranged by the file name order

This script will also create one helpful bash script invoking `text2workspace` program, which can be used on machines with HiggsCombine set up.

Normally, to save time, other frameworks may generate the intermediate 2D histogram templates (containing jet pT versus jet mass distribution) for fast datacard generation in case the user wants to adjust the jet pT range. Unfortunately this may lead to bugs since the 2D histogram may not always have the exact pT ranges encoded. To avoid this surprise, **this script will only generate 1D distribution and no intermediate 2D histogram templates**. 

Furthermore, to offer more flexibility in event categories which may not entirely rely on one pT variable only (such as scale factor measurements for two or more variables, where the categories do not have to follow in the grid fashion), instead of only defining pT ranges, **you can (and must) define your own event categories**. This means, for each event category, you must include all the variables needed in the rule associated with the category.

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
- `categories`: _Tagging categories_ to be added to the datacard. (This is equivalent to _processes_ in HiggsCombine.) Each key represents a tagging category, and must contain the list of processes (under `processes` key) and further cuts for that category (under `cut` key).
- `tagger`: Details on the designated tagger.
    - `name`: Name of the tagger to be used in files
    - `varname`: Name of the discrimnator _as seen in the ntuple files_.
    - `cut`: Discriminator cut for the tagger. Passing and failing events are defined as events with discriminator values higher and lower than this cut respectively.
- `distribution`: Details on the final distribution for scale factor measurement.
    - `mass_variable`: Target variable to be populated in the distribution. Usually jet mass is used.
    - `mass_range`: Range of the target variable.
    - `mass_bins`: Number of bins in the distribution.
    - `event_categories`: Definitions for _event categories_. Every event should be categorised into one of the event categories (provided that they are orthogonal), and then further classified into passing and failing categories based on the tagger. (This is equivalent to _bins_ in HiggsCombine.) **Required.** 
    
      Each category should contain the following keys:
        - `name`: Name of the event category. **The name used here will be used as output file names.**
        - `rule`: Rule of the event category. This rule will be plugged directly into `TTree.Project` method, so the rules defined here should be in the compatible ROOT format.
- `uncertainties`: Details on uncertainties to be added into the datacard. Each key is a uncertainty name, and should contain `mode` key inside. Currently `mode` supports `lnN`, `factor`, and `file`, and has different behaviours as follows:
    - `lnN`: Log-normal uncertainty. If no `category` key is present, this uncertainty will be applied to all _tagging categories_ with the specified `size` value. Specify tagging categories to apply this uncertainty using `category` key.
    - `factor`: Shape uncertainty calculated from _nominal_ input files with the designated expression. Must contain keys `up` and `down`.
    - `file`: Shape uncertainty calculated from files specified in `unc_files`. In this case, the uncertainty name must be the same as specified in `processes` key.
- `perfileweights`: _(Optional)_ Adds a new branch with a _constant_ value to ROOT files. New branches will be added directly to the TTree named in `treename` in specified ROOT files before any histograms are made, which is ideal for updating the cross section for certain files. This option should contain a list of dictionaries following this pattern:
    - `name`: Name for a new branch of the target TTree (named in `treename`)
    - `value`: Value of a new branch. Must be constant number. _Values calculated based on other columns are not supported._
    - `files`: List of ROOT file paths for the new branch to be added.

  See `specfile_test_2022.yaml` for examples on how to use this option.

#### `plot_histograms.py`
The YAML input for `plot_histograms.py` has a different structure, aimed at plotting the prefit and postfit histograms as follows:
- `lumi`: Integrated luminosity. **Required.**
- `xlabel`: Label on the x-axis. **Required.** The range of x-axis is automatically determined by the histogram bins.
- `savedir`: Save directory **Required**
- `categories`: _Tagging categories_ as defined in the YAML input file for `make_histograms.py`. However, each key (name of the tagging category) must contain the following sub-keys:
    - `color`: Histogram colour, used in both prefit and postfit plots
    - `propername`: Name to be shown in the legend
- `eventcats`: _List_ of event categories to be plotted. Must contain the following sub-keys:
    - `name`: Event category name as defined in `event_categories` option in the input YAML file for `make_histograms.py`
    - `propername`: Text to be added in the legend of the plot, _not required_
    - `prefitfile`: Path to prefit file, or the input ROOT file generated from `make_histograms.py`
    - `postfitfile`: Path to postfit file, or the output ROOT file from FitDiagnostics method of Higgs Combine

### What's inside the output ROOT file
The output ROOT file contains _all_ 1D histograms including passing and failing distributions. The naming convention is as follows:
- For data histograms, `data_{EVENT_CATEGORY}_pass` or `data_{EVENT_CATEGORY}_fail`
- For MC histograms, 
    - `{TAGGING_CATEGORY}_{EVENT_CATEGORY}_pass_nominal`
    - `{TAGGING_CATEGORY}_{EVENT_CATEGORY}_pass_{UNCERTAINTY}Up`
    - `{TAGGING_CATEGORY}_{EVENT_CATEGORY}_pass_{UNCERTAINTY}Down`
    - `{TAGGING_CATEGORY}_{EVENT_CATEGORY}_fail_nominal`
    - `{TAGGING_CATEGORY}_{EVENT_CATEGORY}_fail_{UNCERTAINTY}Up`
    - `{TAGGING_CATEGORY}_{EVENT_CATEGORY}_fail_{UNCERTAINTY}Down`

If `--diagnosis` option is turned on for `make_histograms.py`, another ROOT file, per event category, will be created with the name `diagnosis_{EVENT_CATEGORY}.root` The naming convention in this file is 
`{PROCESS}_{UNCERTAINTY}_{up|down}_{FILE_INDEX}_{EVENT_CATEGORY}_{TAGGING_CATEGORY}_{pass|fail}`
where `FILE_INDEX` represents the order of the input file specified in the YAML file.