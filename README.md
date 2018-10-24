# SKFlatMaker

## Environment
```bash
# Setup for >= CMSSW_9_4_9_cand2
export SCRAM_ARCH=slc6_amd64_gcc630

# For machines with CVMFS
export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
source $VO_CMS_SW_DIR/cmsset_default.sh

# Make CMSSW directory
cmsrel CMSSW_9_4_9_cand2
cd CMSSW_9_4_9_cand2/src
cmsenv

# EGamma smearing
# https://twiki.cern.ch/twiki/bin/view/CMS/EgammaMiniAODV2#2017_MiniAOD_V2
git cms-init
git cms-merge-topic cms-egamma:EgammaPostRecoTools_940 #just adds in an extra file to have a setup function to make things easier

# https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2#Recipe_for_regular_users_formats
# The Fall17_94X_V2 ID modules are included by default for CMSSW_10_2_X and higher, for CMSSW_9_4_X you can obtain them through following git cms-merge-topic: 
# Both CutBased and MVA
git cms-merge-topic guitargeek:EgammaID_9_4_X
scram b -j 8

# MET EE Noise filter
# https://hypernews.cern.ch/HyperNews/CMS/get/JetMET/1865.html
# https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETUncertaintyPrescription#Instructions_for_9_4_X_X_9_for_2
git cms-merge-topic cms-met:METFixEE2017_949_v2

# Copy this code
git clone git@github.com:CMSSNU/SKFlatMaker.git
git checkout <branch or tag>

# Compile
scram b -j 8

# Setup
cd $CMSSW_BASE/src/SKFlatMaker
source setup.sh
```

# test runs
```bash
cd $CMSSW_BASE/src/SKFlatMaker/SKFlatMaker/test/
cmsRun RunSKFlatMaker.py year=2016 sampletype=DATA ## Run 2016 DATA
cmsRun RunSKFlatMaker.py year=2016 sampletype=MC ## Run 2016 MC
cmsRun RunSKFlatMaker.py year=2017 sampletype=DATA ## Run 2017 DATA
cmsRun RunSKFlatMaker.py year=2017 sampletype=MC ## Run 2017 MC
```
