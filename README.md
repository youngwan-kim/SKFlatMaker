# SKFlatMaker
* SKFlat Entuple maker for physics analysis with 2017 CMS data
   *  ntuple ()

## Environment
	cmsrel CMSSW_9_4_2
   	cd CMSSW_9_4_2/src
   	cmsenv
   	git cms-init
   	#E/gamma related setting
   	#cut based ID
   	git cms-merge-topic lsoffi:CMSSW_9_4_0_pre3_TnP
   	#MVA based ID
   	git cms-merge-topic guitargeek:ElectronID_MVA2017_940pre3
   	scram b -j 9
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