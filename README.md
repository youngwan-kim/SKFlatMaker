# SKFlatMaker
* SKFlat Entuple maker for physics analysis with 2017 CMS data
   *  ntuple ()


## Environment
```bash
# setup for CMSSW_9_4_X
export SCRAM_ARCH=slc6_amd64_gcc630 #Only for machines with different default compiler. Applying this line at lxplus might cause crash or error during scram

# For machines with CVMFS
export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
source $VO_CMS_SW_DIR/cmsset_default.sh

# make CMSSW directory
cmsrel CMSSW_9_4_4
cd CMSSW_9_4_4/src
cmsenv

# add packages
git cms-init
# based on wiki page https://twiki.cern.ch/twiki/bin/view/CMS/Egamma2017DataRecommendations#How_to_run_the_Scale_Smearing_co
git cms-merge-topic lsoffi:CMSSW_9_4_0_pre3_TnP
git cms-merge-topic guitargeek:ElectronID_MVA2017_940pre3
git cms-merge-topic cms-egamma:EGIDV1AndScaleSmear_940
git cms-merge-topic cms-egamma:EgammaPostRecoTools_940 #just adds in an extra file to have a setup function to make things easier 

scram b -j 8

# Add the area containing the MVA weights (from cms-data, to appear in “external”).
# Note: the “external” area appears after “scram build” is run at least once, as above
#
cd $CMSSW_BASE/external
# below, you may have a different architecture, this is just one example from lxplus
cd slc6_amd64_gcc630/
git clone https://github.com/lsoffi/RecoEgamma-PhotonIdentification.git data/RecoEgamma/PhotonIdentification/data
cd data/RecoEgamma/PhotonIdentification/data
git checkout CMSSW_9_4_0_pre3_TnP
cd $CMSSW_BASE/external
cd slc6_amd64_gcc630/
git clone https://github.com/lsoffi/RecoEgamma-ElectronIdentification.git data/RecoEgamma/ElectronIdentification/data
cd data/RecoEgamma/ElectronIdentification/data
git checkout CMSSW_9_4_0_pre3_TnP
# Go back to the src/
cd $CMSSW_BASE/src

# now we need to get the .dat files for the scale and smearing
cd $CMSSW_BASE/external
# below, you may have a different architecture, this is just one example from lxplus
cd slc6_amd64_gcc630/
git clone git@github.com:Sam-Harper/EgammaAnalysis-ElectronTools.git data/EgammaAnalysis/ElectronTools/data
cd data/EgammaAnalysis/ElectronTools/data
git checkout ReReco17NovScaleAndSmearing 
# Go back to the src/
cd $CMSSW_BASE/src

# compile
scram b -j 8

# copy this code
git clone https://github.com/sungbinoh/SKFlatMaker.git -b <branch_name>

# compile
scram b -j 8

# for crab submission, we need to remove old egamma files
cd $CMSSW_BASE/external
rm -rf slc6_amd64_gcc630/data/RecoEgamma/ElectronIdentification/data/Spring1*
rm -rf slc6_amd64_gcc630/data/RecoEgamma/ElectronIdentification/data/PHYS14
rm -rf slc6_amd64_gcc630/data/RecoEgamma/ElectronIdentification/data/*.xml

rm -rf slc6_amd64_gcc630/data/RecoEgamma/PhotonIdentification/data/Spring1*
rm -rf slc6_amd64_gcc630/data/RecoEgamma/PhotonIdentification/data/PHYS14



# test cmsRun
cd SKFlatMaker/SKFlatMaker/ntuples/suoh_test
# modify DATA_cfg_test_2017promptReco.py file, eg) TESTFILE_DATA for your test rootfile, isMC also
voms-proxy-init --voms cms
cmsRun DATA_test_Nov17_ReReco.py #or DATA_test_2017_PromptReco.py, MC_test_944_mc2017.py
```
