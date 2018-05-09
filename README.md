# SKFlatMaker
* SKFlat Entuple maker for physics analysis with 2017 CMS data
   *  ntuple ()


## Environment
```bash
# setup for >= 9_4_6_patch1
export SCRAM_ARCH=slc6_amd64_gcc630 #Only for machines with different default compiler. Applying this line at lxplus might cause crash or error during scram

# For machines with CVMFS
export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
source $VO_CMS_SW_DIR/cmsset_default.sh

# make CMSSW directory
cmsrel 9_4_6_patch1
cd 9_4_6_patch1/src
cmsenv

# add packages
git cms-init
# based on wiki page https://twiki.cern.ch/twiki/bin/view/CMS/Egamma2017DataRecommendations#How_to_run_the_Scale_Smearing_co
git cms-merge-topic cms-egamma:EgammaPostRecoTools_940 #just adds in an extra file to have a setup function to make things easier

# copy this code
git clone https://github.com/sungbinoh/SKFlatMaker.git -b <branch_name>

# compile
scram b -j 8

# test cmsRun
cd SKFlatMaker/SKFlatMaker/ntuples/suoh_test
# modify DATA_cfg_test_2017promptReco.py file, eg) TESTFILE_DATA for your test rootfile, isMC also
voms-proxy-init --voms cms
cmsRun DATA_test_Nov17_ReReco.py #or DATA_test_2017_PromptReco.py, MC_test_944_mc2017.py
```
