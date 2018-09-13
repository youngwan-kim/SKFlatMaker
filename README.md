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
scram b -j 8

# Copy this code
git clone git@github.com:jedori0228/SKFlatMaker.git
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
