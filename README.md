# SKFlatMaker
## Release notes
### Run2UltraLegacy_v3 (2022 June)
* Move to UL MiniAOD v2. Add tau, puppiMET, etc. See more details [here](https://github.com/CMSSNU/SKFlatMaker/issues/63)
* Processed samples are listed in below text files. But please also check the HN google sheet for more updated list  
  * [2016preVFP_DATA](SKFlatMaker/script/CRAB3/2016preVFP_DATA.txt), [2016preVFP_MC](SKFlatMaker/script/CRAB3/2016preVFP_MC.txt)  
  * [2016postVFP_DATA](SKFlatMaker/script/CRAB3/2016postVFP_DATA.txt), [2016postVFP_MC](SKFlatMaker/script/CRAB3/2016postVFP_MC.txt)  
  * [2017_DATA](SKFlatMaker/script/CRAB3/2017_DATA.txt), [2017_MC](SKFlatMaker/script/CRAB3/2017_MC.txt)  
  * [2018_DATA](SKFlatMaker/script/CRAB3/2018_DATA.txt), [2018_MC](SKFlatMaker/script/CRAB3/2018_MC.txt)  
  
## Setup environment
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
git checkout -b Run2UltraLegacy_v3 Run2UltraLegacy_v3 #### use the tag

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

<<<<<<< HEAD
## Submitting crab jobs
1. Make a sample list to be processed at *SKFlatMaker/script/CRAB3* directory. **The filename should include the era.** See [SKFlatMaker/script/CRAB3/2018_DATA.txt](SKFlatMaker/script/CRAB3/2018_DATA.txt) for example.
2. (optional) The number of events in the second column are not necessary but *checkDAS.py* script will can help to get the numbers. If a sample is VALID state, it will print the number of events. If a sample is PRODUCTION state, it will print the progress. It will overwrite the sample list file with `-w` option.
```
python checkDAS.py 2018_MYSAMPLELIST.txt
python checkDAS.py 2018_MYSAMPLELIST.txt -w
```
3. (optional) If you want to save events weight for the scale variation (LHE level), PDF reweight (LHE level), parton shower weight (Pythia level) etc, You should provide the weight IDs you want to store (weight map file). First, get logs containing weight description using `getLog.py` script and then make weight map files using `parseLog.py` script. See examples in the relevant directories and [files](SKFlatMaker/script/Weight/data/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8__RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1.txt).
```
cd $CMSSW_BASE/src/SKFlatMaker/SKFlatMaker/script/Weight
python getLog.py ../CRAB3/2018_MYSAMPLELIST.txt
python parseLog.py logs/
```
4. Make crab configuration files using `MakeCrab.py` script. There will be warnings if you don't provide a weight map file or it miss some basic weights such as Scale, PDF, AlphaS weight. In that case these weights will not be stored.
```
cd $CMSSW_BASE/src/SKFlatMaker/SKFlatMaker/script/CRAB3
python MakeCrab.py 2018_MYSAMPLELIST.txt
```
5. Now submit crab jobs
```
cd $SKFlatTag/2018/crab_submission_MC/
##submit specific sample
crab submit -c SubmitCrab__MYSAMPLE.py
##submit all
ls SubmitCrab*.py|xargs -i crab submit -c {} 
```

## Adding a variable
1. Define the variable in *SKFlatMaker/interface/SKFlatMaker.h*
2. Initialize or clear the variable at *SKFlatMaker::analyze* function in *SKFlatMaker/src/SKFlatMaker.cc*
```c++
//SKFlatMaker::analyze
myvar = -1; // for primitive type
myvector.clear(); // for vector type
```
3. Add branch to the tree
```c++
//SKFlatMaker::beginJob
DYTree->Branch("myvar",&myvar,"myvar/F"); // for primitive type
DYTree->Branch("myvector", "vector<float>", &myvector); // for vector type
```
4. Fill the value at proper *SKFlatMaker::fill??* function
```c++
//SKFlatMaker::fillSOMETHING
myvar=VALUE; // for primitive type
myvector.push_back(VALUE); // for vector type
```
5. Insert the variable on proper varlist file at *SKFlatMaker/script/VarList/*. The first column is type, the second is variable name and the last is 0 for primitive type and 1 for vector type.
6. Check using *check_varlist.py* script
```
$CMSSW_BASE/src/SKFlatMaker/SKFlatMaker/script/VarList
python check_varlist.py
```
