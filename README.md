# TopSF
A _better_ framework for scale factor measurement, designed with flexibility in mind.

**[We have a Wiki now! Yay!](../../wiki)**

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