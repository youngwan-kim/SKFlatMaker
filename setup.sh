#!/bin/bash
if [ -z "$CMSSW_BASE" ]; then
    if [ -e "../../src" ]; then
	source /cvmfs/cms.cern.ch/cmsset_default.sh
	cd ../../src
	eval `scramv1 runtime -sh`
	cd -
    else
	echo "Wrong directory structure" >&2
	exit 1
    fi
fi
export SKFlatTag=Run2UltraLegacy_v3
export SKFlatWD=$CMSSW_BASE/src/SKFlatMaker/
