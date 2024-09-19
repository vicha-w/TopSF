# How to use the generated HiggsCombine datacard and accompanying ROOT file
`make_histograms.py` accepts the YAML specfile indicating the locations of input ROOT files and all configurations, and outputs one HiggsCombine datacard and one ROOT file (containing required histograms) per pT category. Here are the instructions on how to use this datacard for statistical measurements using HiggsCombine.

1. **Install `TagAndProbeExtended` model into HiggsCombine.** Scale factor measurements require the tag-and-probe model, which is described in `TagAndProbeExtended` model. The model is by default not available, so you need to copy and compile it.
    1. Download the model file from https://github.com/cms-jet/ParticleNetSF/blob/ParticleNet_TopW_SFs_NanoV9/TagAndProbeExtended.py
    2. Assuming that you have already installed HiggsCombine, copy the model file to `HiggsAnalysis/CombinedLimit/python`
    3. **Compile** the model file with `scram build`
2. **Convert the datacard to a workspace file.** For portability, converting the datacard to a workspace using `text2workspace.py` is strongly recommended. `make_histograms.py` also generates the `text2workspace.py` script required for this conversion, aptly named `combine_script.sh`. The output files for this conversion are `workspace_*.root`, as shown in the `text2workspace.py` command:
    ```bash
    text2workspace.py -m 125 -P HiggsAnalyis.CombinedLimit.TagAndProbeExtended:tagAndProbe <PT_RANGE>.txt -o workspace_<PT_RANGE>.root --PO=categories=<TAGGING_CATEGORIES>
    ```
    You can change the workspace name after the option `-o`. The resulting workspace files will have as many scale factor parameters as the number of categories specified in the option `--PO`, and will have `SF_` preceding the names of categories. For example, if you have four categories named `topmatched`, `wmatched`, `nonmatched`, and `other`, the scale factors will be named `SF_topmatched`, `SF_wmatched`, `SF_nonmatched`, and `SF_other` respectively.

    The tag-and-probe model we are using will define four scale factors with the range of [0, 2]. See [Tag-and-probe model used](#tag-and-probe-model-used) section for more details.

    **All `SF_*` parameters are declared as parameters of interest by this model by default.** If you want to focus on only one parameter, such as top-tagging scale factor `SF_topmatched`, add the option `--redefineSignalPOIs` into `combine`, e.g. `--redefineSignalPOIs SF_topmatched`.
3. **Use the workspace file as HiggsCombine input.** Once the workspace file is created, you now have a creative freedom on whatever you want to measure with it. Refer to the example combine script for ideas to get you started.

## MC/Data normalisation factor
For the scale factor measurement, you may want to ensure that the total yield of MC events matches the total yield of data events in one fit. This is done using `norm_match_mc_data` parameter, defined at the bottom of the datacard just before `autoMCstats`. This parameter is calculated automatically from data yield divided by MC yield, and are frozen by default with `nuisance edit freeze norm_match_mc_data`. By converting the datacard into a workspace file, this parameter is also included in the workspace file.

## HiggsCombine example commands
In the following code, assume that `workspace.root` file is our workspace file. Extra fit options such as `--cminDefaultMinimizerTolerance` can be added as indicated. Refer to HiggsCombine documentation for more details.
- Measuring top-tagging scale factor. **Note that output from FitDiagnostics is required for postfit plots.**
    ```bash
    # Just the scale factor, please
    combine -M MultiDimFit workspace.root --redefineSignalPOIs SF_topmatched <EXTRA FIT OPTIONS HERE>

    # with FitDiagnostics output file, please
    combine -M FitDiagnostics workspace.root --redefineSignalPOIs SF_topmatched <EXTRA FIT OPTIONS HERE>
    ```
- One dimensional NLL scan
    ```bash
    OUTPUTNAME=topsf_scan_nll

    combine -M MultiDimFit workspace.root --algo=grid --points 201 --alignEdges 1 --setParameterRanges SF_topmatched=0,2 --redefineSignalPOIs SF_topmatched --name ${OUTPUTNAME} <EXTRA FIT OPTIONS HERE>

    plot1DScan.py higgsCombine${OUTPUTNAME}.MultiDimFit.mH120.root --POI SF_topmatched
    ```
- Impacts plot
    ```bash
    # Initial fit
    combineTool.py -M Impacts -d workspace.root -m 125 --doInitialFit --redefineSignalPOIs SF_topmatched <EXTRA FIT OPTIONS HERE>

    # Fits for all parameters
    # --exclude 'rgx{prop.*}' excludes all statistical uncertainties
    combineTool.py -M Impacts -d workspace.root -m 125 --doFits --exclude 'rgx{prop.*}' --redefineSignalPOIs SF_topmatched <EXTRA FIT OPTIONS HERE>

    # Summarise the results into JSON file
    combineTool.py -M Impacts -d workspace.root -m 125 -o impacts.json --exclude 'rgx{prop.*}' --redefineSignalPOIs SF_topmatched
    ```

## Tag-and-probe model used
The tag-and-probe model used in this framework is available at https://github.com/cms-jet/ParticleNetSF/blob/ParticleNet_TopW_SFs_NanoV9/TagAndProbeExtended.py.

In the default physics model used in HiggsCombine, there will be one parameter of interest, which is the signal strength $\\mu$ (or `r` in HiggsCombine fitting results). Signal strength is the amount of how much the signal process event yield actually appears in the data distribution compared to MC prediction, and scales the signal yield across all event categories by its value. This model, however, is not suitable for scale factor measurements.

In this custom tag-and-probe model, each _tagging category_ will have a designated scale factor with the prefit value of 1 and a range of [0, 2]. The yield for each tagging category is scaled differently for passing and failing distributions as follows:
- $P \\rightarrow P' = P \\times SF$ 
- $F \\rightarrow F' = P + F - P \\times SF$

where $P$ and $F$ represents the passing and failing event yields, and $SF$ is the scale factor. This custom scaling ensures that the total event yield between passing and failing events is preserved. The model, however, also ensures that the failing event yield must never be below zero."