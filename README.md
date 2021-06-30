# SKFlatMaker

## Environment
```bash
export SCRAM_ARCH=slc7_amd64_gcc700

# For machines with CVMFS
export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
source $VO_CMS_SW_DIR/cmsset_default.sh

# Make CMSSW directory

#### 1) For a test job or developement
cmsrel CMSSW_10_6_26
cd CMSSW_10_6_26/src

#### 2) For production, let's not use the working directory but use new and clean directory
#### Also, I recommend using lxplus
scram p -n Run2UltraLegacy_v1__CMSSW_10_6_26 CMSSW CMSSW_10_6_26
cd Run2UltraLegacy_v1__CMSSW_10_6_26/src

#### Then,
cmsenv
git cms-init

######################
#### EGamma smearing
#### https://twiki.cern.ch/twiki/bin/view/CMS/EgammaUL2016To2018
######################

git cms-addpkg RecoEgamma/EgammaTools  ### essentially just checkout the package from CMSSW
git clone https://github.com/cms-egamma/EgammaPostRecoTools.git
mv EgammaPostRecoTools/python/EgammaPostRecoTools.py RecoEgamma/EgammaTools/python/.
git clone https://github.com/jainshilpi/EgammaAnalysis-ElectronTools.git -b UL2018 EgammaAnalysis/ElectronTools/data/
git cms-addpkg EgammaAnalysis/ElectronTools
scram b -j 4

######################
#### Now SKFlatMaker
######################

# Copy this code
git clone https://github.com/CMSSNU/SKFlatMaker.git
#Or git clone git@github.com:CMSSNU/SKFlatMaker.git
cd SKFlatMaker

#### 1) For a test job or development 
git checkout master
git checkout -b <testbranch>

#### 2) For production
git checkout Run2UltraLegacy_v1 #### use the tag

cd $CMSSW_BASE/src

# Compile
scram b -j 4

# Setup
cd $CMSSW_BASE/src/SKFlatMaker
source setup.sh

# Now, submitting crab jobs.
cd $CMSSW_BASE/src/SKFlatMaker/SKFlatMaker/script/CRAB3
# edit MakeCrab.py
# txtfilename = "<txtfile which has list of samples you want to run>"
# e.g., if txtfilename = "2018_DATA.txt"
python MakeCrab.py
# then it will print crab submission commands
# copy them somewhere
cd $SKFlatTag/2018/crab_submission_DATA/
# now, run the submission commands
```

# test runs
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
