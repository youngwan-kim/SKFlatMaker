# SKFlatMaker

## Environment
```bash
export SCRAM_ARCH=slc7_amd64_gcc700
source /cvmfs/cms.cern.ch/cmsset_default.sh

# Make CMSSW directory
scram p -n Run2UltraLegacy_v3__CMSSW_10_6_29 CMSSW CMSSW_10_6_29
cd Run2UltraLegacy_v3__CMSSW_10_6_29/src
cmsenv
git cms-init

######################
#### EGamma smearing
#### https://twiki.cern.ch/twiki/bin/view/CMS/EgammaUL2016To2018
######################

git cms-addpkg RecoEgamma/EgammaTools  ### essentially just checkout the package from CMSSW
git clone https://github.com/cms-egamma/EgammaPostRecoTools.git
mv EgammaPostRecoTools/python/EgammaPostRecoTools.py RecoEgamma/EgammaTools/python/.
git clone -b ULSSfiles_correctScaleSysMC https://github.com/jainshilpi/EgammaAnalysis-ElectronTools.git EgammaAnalysis/ElectronTools/data/
git cms-addpkg EgammaAnalysis/ElectronTools
scram b -j 4

######################
#### Now SKFlatMaker
######################

# Copy this code
git clone git@github.com:CMSSNU/SKFlatMaker.git
cd SKFlatMaker
git checkout Run2UltraLegacy
#git checkout -b Run2UltraLegacy_v3 Run2UltraLegacy_v3 #### use the tag


# Compile
cd $CMSSW_BASE/src
scram b -j 4

# Setup
cd $CMSSW_BASE/src/SKFlatMaker
source setup.sh
```

## Test runs
```bash
cd $CMSSW_BASE/src/SKFlatMaker/SKFlatMaker/test/
cmsRun RunSKFlatMaker.py era=2016preVFP sampletype=DATA maxEvents=1000  ## Run 2016a DATA
cmsRun RunSKFlatMaker.py era=2016preVFP sampletype=MC maxEvents=1000    ## Run 2016a MC
cmsRun RunSKFlatMaker.py era=2016postVFP sampletype=DATA maxEvents=1000 ## Run 2016b DATA
cmsRun RunSKFlatMaker.py era=2016postVFP sampletype=MC maxEvents=1000   ## Run 2016b MC
cmsRun RunSKFlatMaker.py era=2017 sampletype=DATA maxEvents=1000        ## Run 2017 DATA
cmsRun RunSKFlatMaker.py era=2017 sampletype=MC maxEvents=1000          ## Run 2017 MC
cmsRun RunSKFlatMaker.py era=2018 sampletype=DATA maxEvents=1000        ## Run 2018 DATA
cmsRun RunSKFlatMaker.py era=2018 sampletype=MC maxEvents=1000          ## Run 2018 MC
```
Some useful options  
- weightmap: weightmap file path or string: ex) Scale[1,9],PDF[1001,1100],AlphaS[1101,1102],AlphaSScale[1.5,1.5]  
  weightmap files are stored at `$SKFlatWD/SKFlatMaker/script/Weight/data`  
- debug: Set print level defalt=0  
- `python RunSKFlatMaker.py help` for other options  

## Sample processing
1. make a sample list you want to process at `$CMSSW_BASE/src/SKFlatMaker/SKFlatMaker/script/CRAB3`. The file name should include the correct era.
2. (optional) check the sample status using `checkDAS.py` script. If a sample is VALID state, it will print the number of events. If a sample is PRODUCTION state, it will print the progress. It will overwrite the sample list file with `-w` option.
```
python checkDAS.py SAMPLELISTFILE
python checkDAS.py SAMPLELISTFILE -w
```
3. (optional) If you want to save events weight for the scale variation (LHE level), PDF reweight (LHE level), parton shower weight (Pythia level) etc, You should provide the weight IDs you want to store (weight map file). First, get logs containing weight description using `getLog.py` script and then make weight map files using `parseLog.py` script. See examples in the relevant directories.
```
cd $CMSSW_BASE/src/SKFlatMaker/SKFlatMaker/script/Weight
python getLog.py ../CRAB3/SAMPLEFILELIST
python parseLog.py logs/
```
4. Make crab script using `MakeCrab.py` script. There will be warnings if you don't provide a weight map file or it miss some basic weights such as Scale, PDF, AlphaS weight. In that case these weights will not be stored.
```
cd $CMSSW_BASE/src/SKFlatMaker/SKFlatMaker/script/CRAB3
python MakeCrab.py SAMPLELISTFILE
```
5. Now, submitting crab jobs. For example, 2016preVFP DY MC
```
cd $CMSSW_BASE/src/SKFlatMaker/SKFlatMaker/script/CRAB3/$SKFlatTag/2016preVFP/crab_submission_MC/
crab submit -c SubmitCrab__DYJetsToEE_M-50_massWgtFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos__RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1.py
```

## Processed samples
* [2016preVFP_DATA](SKFlatMaker/script/CRAB3/2016preVFP_DATA.txt)
* [2016preVFP_MC](SKFlatMaker/script/CRAB3/2016preVFP_MC.txt)

## Release note
v3 (2022Jun)
* MiniAODv1 -> MiniAODv2
* Add Tau object, Puppi MET, and other variables. See https://github.com/CMSSNU/SKFlatMaker/issues/63 for the details
