# SKFlatMaker

## Environment
```bash
# Setup for >= CMSSW_10_2_5
export SCRAM_ARCH=slc7_amd64_gcc700

# For machines with CVMFS
export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
source $VO_CMS_SW_DIR/cmsset_default.sh

# Make CMSSW directory

#### 1) For a test job or developement
cmsrel CMSSW_10_4_0_patch1
cd CMSSW_10_4_0_patch1/src

#### 2) For production, let's not use the working directory but use new and clean directory
#### Also, I recommend using lxplus
scram p -n Run2Legacy_v3__CMSSW_10_4_0_patch1 CMSSW CMSSW_10_4_0_patch1
cd Run2Legacy_v3__CMSSW_10_4_0_patch1/src

#### Then,
cmsenv
git cms-init

######################
#### EGamma smearing
#### https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaPostRecoRecipes#2018_Data_MC
######################

git cms-merge-topic cms-egamma:EgammaPostRecoTools #just adds in an extra file to have a setup function to make things easier 

#### Below gives conflict
#git cms-merge-topic cms-egamma:PhotonIDValueMapSpeedup1029 #optional but speeds up the photon ID value module so things fun faster
#### Below gives conflict
#git cms-merge-topic cms-egamma:slava77-btvDictFix_10210 #fixes the Run2018D dictionary issue, see https://github.com/cms-sw/cmssw/issues/26182, may not be necessary for later releases, try it first and see if it works

#now to add the scale and smearing for 2018 (eventually this will not be necessary in later releases but is harmless to do regardless)
git cms-addpkg EgammaAnalysis/ElectronTools
rm EgammaAnalysis/ElectronTools/data -rf
git clone git@github.com:cms-data/EgammaAnalysis-ElectronTools.git EgammaAnalysis/ElectronTools/data
#now build everything
scram b -j 8

###############################
#### Rerun ecalBadCalibfilter
#### https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2#How_to_run_ecal_BadCalibReducedM
###############################

git cms-addpkg RecoMET/METFilters

#########################
#### L1Prefire reweight
#### https://twiki.cern.ch/twiki/bin/view/CMS/L1ECALPrefiringWeightRecipe
#########################

git clone git@github.com:jedori0228/L1Prefiring.git
#### Do not compile here, we need SKFlatMaker cloned

######################
#### Now SKFlatMaker
######################

# Copy this code
git clone git@github.com:CMSSNU/SKFlatMaker.git
cd SKFlatMaker

#### 1) For a test job or development 
git checkout master
git checkout -b <testbranch>

#### 2) For production
git checkout Run2Legacy_v3 #### use the tag

cd $CMSSW_BASE/src

# Compile
scram b -j 8

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
cmsRun RunSKFlatMaker.py year=2016 sampletype=DATA maxEvents=1000 ## Run 2016 DATA
cmsRun RunSKFlatMaker.py year=2016 sampletype=MC maxEvents=1000 ## Run 2016 MC
cmsRun RunSKFlatMaker.py year=2017 sampletype=DATA maxEvents=1000 ## Run 2017 DATA
cmsRun RunSKFlatMaker.py year=2017 sampletype=MC maxEvents=1000 ## Run 2017 MC
cmsRun RunSKFlatMaker.py year=2018 sampletype=DATA maxEvents=1000 ## Run 2018 DATA
cmsRun RunSKFlatMaker.py year=2018 sampletype=MC maxEvents=1000 ## Run 2018 MC
```
