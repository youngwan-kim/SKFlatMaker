# SKFlatMaker
* SKFlat Entuple maker for physics analysis with 2017 CMS data
   *  ntuple ()


## Environment
	#setup for CMSSW_9_4_X
	export SCRAM_ARCH=slc6_amd64_gcc630 #Only for machines with different default compiler. Applying this line at lxplus might cause crash or error during scram
	cmsrel CMSSW_9_4_2
   	cd CMSSW_9_4_2/src
   	cmsenv	
				
	#add packages
   	git cms-init
   	#E/gamma related setting(https://twiki.cern.ch/twiki/bin/view/CMS/EgammaIDRecipesRun2)
   	git cms-merge-topic lsoffi:CMSSW_9_4_0_pre3_TnP
   	git cms-merge-topic guitargeek:ElectronID_MVA2017_940pre3
	scram b -j 8 # This will make additional directories 
	cd $CMSSW_BASE/external
   	cd slc6_amd64_gcc630/	
   	git clone https://github.com/lsoffi/RecoEgamma-PhotonIdentification.git data/RecoEgamma/PhotonIdentification/data
   	cd data/RecoEgamma/PhotonIdentification/data
   	git checkout CMSSW_9_4_0_pre3_TnP
   	cd $CMSSW_BASE/external
   	cd slc6_amd64_gcc630/
   	git clone https://github.com/lsoffi/RecoEgamma-ElectronIdentification.git data/RecoEgamma/ElectronIdentification/data
   	cd data/RecoEgamma/ElectronIdentification/data
   	git checkout CMSSW_9_4_0_pre3_TnP
   	cd $CMSSW_BASE/src

	#E/gamma smearing (https://twiki.cern.ch/twiki/bin/view/CMS/EGMSmearer)
	git cms-merge-topic cms-egamma:EGM_94X_v1
	cd EgammaAnalysis/ElectronTools/data
	# download the txt files with the corrections
	git clone https://github.com/ECALELFS/ScalesSmearings.git
	cd ScalesSmearings/
	git checkout Run2017_17Nov2017_v1
	#compile
	cd $CMSSW_BASE/src

	#copy this code
	git clone https://github.com/sungbinoh/SKFlatMaker.git Phys -b MoveToSKFlat

	#compile
	scram b -j 8

	#test cmsRun
	cd Phys/DYntupleMaker/ntuples/suoh_test
	#modify DATA_cfg_test_2017promptReco.py file, eg) TESTFILE_DATA for your test rootfile, isMC also
	voms-proxy-init --voms cms
	cmsRun DATA_test_Nov17_ReReco.py #or DATA_test_2017_PromptReco.py, MC_test_94X_mc2017.py