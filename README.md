# SKFlatMaker

## For 2016 and 2017

### Environment
```bash
# Setup for >= CMSSW_9_4_13
export SCRAM_ARCH=slc6_amd64_gcc630

# For machines with CVMFS
export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
source $VO_CMS_SW_DIR/cmsset_default.sh

# Make CMSSW directory
cmsrel CMSSW_9_4_13
cd CMSSW_9_4_13/src
cmsenv
git cms-init

######################
#### EGamma smearing
#### https://twiki.cern.ch/twiki/bin/view/CMS/EgammaMiniAODV2#2017_MiniAOD_V2
######################

git cms-merge-topic cms-egamma:EgammaPostRecoTools

##########################
#### MET EE Noise filter
#### https://hypernews.cern.ch/HyperNews/CMS/get/JetMET/1865.html
#### https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETUncertaintyPrescription#Instructions_for_9_4_X_X_9_for_2
##########################

git cms-merge-topic cms-met:METFixEE2017_949_v2

###############################
#### Rerun ecalBadCalibfilter
#### https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2#How_to_run_ecal_BadCalibReducedM
###############################

git cms-addpkg RecoMET/METFilters

#########################
#### L1Prefire reweight
#### https://twiki.cern.ch/twiki/bin/view/CMS/L1ECALPrefiringWeightRecipe
#########################

#git cms-merge-topic lathomas:L1Prefiring_9_4_9
git cms-merge-topic jedori0228:L1Prefiring_9_4_9__UseFileFileInPath

################################################################################################
################################################################################################

# Copy this code
git clone git@github.com:CMSSNU/SKFlatMaker.git
cd SKFlatMaker
git checkout <branch or tag>
cd $CMSSW_BASE/src

# Compile
scram b -j 8

# Setup
cd $CMSSW_BASE/src/SKFlatMaker
source setup.sh
```
### For the Production
I recommend you to make a clean working directory when we do the central production.

Also, I recommend using lxplus

```bash
# scram p -n <directory name> CMSSW <cmssw release>
# below, assuming we are using tag:Run2Legacy_v3, in cmssw:CMSSW_9_4_13
export SCRAM_ARCH=slc6_amd64_gcc630
scram p -n Run2Legacy_v3__CMSSW_9_4_13 CMSSW CMSSW_9_4_13
cd Run2Legacy_v3__CMSSW_9_4_13/src
cmsenv
# crab setup
source /cvmfs/cms.cern.ch/crab3/crab.sh

###############################################
# Environmenet setup. Basically same as above
###############################################

git cms-init

######################
#### EGamma smearing
#### https://twiki.cern.ch/twiki/bin/view/CMS/EgammaMiniAODV2#2017_MiniAOD_V2
######################

git cms-merge-topic cms-egamma:EgammaPostRecoTools

##########################
#### MET EE Noise filter
#### https://hypernews.cern.ch/HyperNews/CMS/get/JetMET/1865.html
#### https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETUncertaintyPrescription#Instructions_for_9_4_X_X_9_for_2
##########################

git cms-merge-topic cms-met:METFixEE2017_949_v2

###############################
#### Rerun ecalBadCalibfilter
#### https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2#How_to_run_ecal_BadCalibReducedM
###############################

git cms-addpkg RecoMET/METFilters

#########################
#### L1Prefire reweight
#### https://twiki.cern.ch/twiki/bin/view/CMS/L1ECALPrefiringWeightRecipe
#########################

#git cms-merge-topic lathomas:L1Prefiring_9_4_9
git cms-merge-topic jedori0228:L1Prefiring_9_4_9__UseFileFileInPath

# Now SKFlatMaker

git clone git@github.com:CMSSNU/SKFlatMaker.git
cd SKFlatMaker
git checkout Run2Legacy_v3 ## This let us use exactly same tag without any modifications..
source setup.sh
cd $CMSSW_BASE/src
scram b -j 8

# Now, submitting crab jobs.
cd $CMSSW_BASE/src/SKFlatMaker/SKFlatMaker/script/CRAB3
# txtfilename = "<txtfile which has list of samples you want to run>"
# e.g., if txtfilename = "2017_DATA.txt"
python MakeCrab.py <txtfilename>
# then it will print crab submission commands
# copy them somewhere
cd $SKFlatTag/2017/crab_submission_DATA/
# now, run the submission commands
```
## For 2018

### Environment
```bash
# Setup for >= CMSSW_10_2_5
export SCRAM_ARCH=slc6_amd64_gcc700

# For machines with CVMFS
export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
source $VO_CMS_SW_DIR/cmsset_default.sh

# Make CMSSW directory
cmsrel CMSSW_10_4_0_patch1
cd CMSSW_10_4_0_patch1/src
cmsenv
git cms-init

######################
#### EGamma smearing
#### https://twiki.cern.ch/twiki/bin/view/CMS/EgammaMiniAODV2#2018_MiniAOD
#### https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaPostRecoRecipes#2018_Preliminary_Energy_Correcti
######################

git cms-merge-topic cms-egamma:EgammaPostRecoTools #just adds in an extra file to have a setup function to make things easier
scram b -j 8

git clone git@github.com:cms-egamma/EgammaAnalysis-ElectronTools.git EgammaAnalysis/ElectronTools/data
cd EgammaAnalysis/ElectronTools/data
git checkout ScalesSmearing2018_Dev
cd -
git cms-merge-topic cms-egamma:EgammaPostRecoTools_dev
scram b -j 8

###############################
#### Rerun ecalBadCalibfilter
#### https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2#How_to_run_ecal_BadCalibReducedM
###############################

git cms-addpkg RecoMET/METFilters

################################################################################################
################################################################################################

# Copy this code
git clone git@github.com:CMSSNU/SKFlatMaker.git
cd SKFlatMaker
git checkout <branch or tag>
cd $CMSSW_BASE/src

# Compile
scram b -j 8

# Setup
cd $CMSSW_BASE/src/SKFlatMaker
source setup.sh
```

### For the Production
I recommend you to make a clean working directory when we do the central production.

Also, I recommend using lxplus

```bash
# scram p -n <directory name> CMSSW <cmssw release>
# below, assuming we are using tag:Run2Legacy_v3, in cmssw:CMSSW_10_4_0_patch1
export SCRAM_ARCH=slc6_amd64_gcc700
scram p -n Run2Legacy_v3__CMSSW_10_4_0_patch1 CMSSW CMSSW_10_4_0_patch1
cd Run2Legacy_v3__CMSSW_10_4_0_patch1/src
cmsenv
# crab setup
source /cvmfs/cms.cern.ch/crab3/crab.sh

###############################################
# Environmenet setup. Basically same as above
###############################################

git cms-init

######################
#### EGamma smearing
#### https://twiki.cern.ch/twiki/bin/view/CMS/EgammaMiniAODV2#2018_MiniAOD
#### https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaPostRecoRecipes#2018_Preliminary_Energy_Correcti
######################

git cms-merge-topic cms-egamma:EgammaPostRecoTools #just adds in an extra file to have a setup function to make things easier
scram b -j 8

git clone git@github.com:cms-egamma/EgammaAnalysis-ElectronTools.git EgammaAnalysis/ElectronTools/data
cd EgammaAnalysis/ElectronTools/data
git checkout ScalesSmearing2018_Dev
cd -
git cms-merge-topic cms-egamma:EgammaPostRecoTools_dev
scram b -j 8

###############################
#### Rerun ecalBadCalibfilter
#### https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2#How_to_run_ecal_BadCalibReducedM
###############################

git cms-addpkg RecoMET/METFilters

# Now SKFlatMaker

git clone git@github.com:CMSSNU/SKFlatMaker.git
cd SKFlatMaker
git checkout Run2Legacy_v3 ## This let us use exactly same tag without any modifications..
source setup.sh
cd $CMSSW_BASE/src
scram b -j 8

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

## test runs
```bash
cd $CMSSW_BASE/src/SKFlatMaker/SKFlatMaker/test/
cmsRun RunSKFlatMaker.py year=2016 sampletype=DATA maxEvents=1000 ## Run 2016 DATA
cmsRun RunSKFlatMaker.py year=2016 sampletype=MC maxEvents=1000 ## Run 2016 MC
cmsRun RunSKFlatMaker.py year=2017 sampletype=DATA maxEvents=1000 ## Run 2017 DATA
cmsRun RunSKFlatMaker.py year=2017 sampletype=MC maxEvents=1000 ## Run 2017 MC
cmsRun RunSKFlatMaker.py year=2018 sampletype=DATA maxEvents=1000 ## Run 2018 DATA
cmsRun RunSKFlatMaker.py year=2018 sampletype=MC maxEvents=1000 ## Run 2018 MC
```
