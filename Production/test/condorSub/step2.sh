#!/bin/bash

# check for incorrect pilot cert
vomsident=$(voms-proxy-info -identity)
echo $vomsident
if [[ $vomsident = *"cmsgli"* ]]; then
	# this is the exit code for "User is not authorized to write to destination site."
	exit 60322
fi

export JOBNAME=""
export PROCESS=""
export OUTDIR=""
export REDIR=""
export OPTIND=1
while [[ $OPTIND -lt $# ]]; do
	# getopts in silent mode, don't exit on errors
	getopts ":j:p:o:x:" opt || status=$?
	case "$opt" in
		j) export JOBNAME=$OPTARG
		;;
		p) export PROCESS=$OPTARG
		;;
		o) export OUTDIR=$OPTARG
		;;
		x) export REDIR=$OPTARG
		;;
		# keep going if getopts had an error
		\? | :) OPTIND=$((OPTIND+1))
		;;
	esac
done

echo "parameter set:"
echo "OUTDIR:     $OUTDIR"
echo "JOBNAME:    $JOBNAME"
echo "PROCESS:    $PROCESS"
echo "REDIR:      $REDIR"
echo ""

# link files from CMSSW dir
ln -s ${CMSSWVER}/src/TreeMaker/Production/test/data
ln -s ${CMSSWVER}/src/TreeMaker/Production/test/runMakeTreeFromMiniAOD_cfg.py

# run CMSSW
ARGS=$(cat args_${JOBNAME}_${PROCESS}.txt)
if [[ -n "$REDIR" ]]; then
	ARGS="$ARGS redir=${REDIR}"
fi
THREADS=$(getFromClassAd RequestCpus)
if [[ -n "$THREADS" ]]; then
	ARGS="$ARGS threads=${THREADS}"
fi
echo "cmsRun runMakeTreeFromMiniAOD_cfg.py ${ARGS} 2>&1"
cmsRun runMakeTreeFromMiniAOD_cfg.py ${ARGS} 2>&1

CMSEXIT=$?

rm runMakeTreeFromMiniAOD_cfg.py

if [[ $CMSEXIT -eq 77 ]]; then
  echo "file already processed"
  exit 0
fi

if [[ $CMSEXIT -ne 0 ]]; then
  rm *.root
  echo "exit code $CMSEXIT, skipping xrdcp"
  exit $CMSEXIT
fi

# check for incorrect pilot cert
vomsident = $(voms-proxy-info -identity)
echo $vomsident
if [[ $vomsident = *"cmsgli"* ]]; then
	# this is the exit code for "User is not authorized to write to destination site."
	rm *.root
	echo "exit code 60322, skipping xrdcp"	
	exit 60322
fi

# write short test script to check if output file has track collection:
echo "
import os, sys, glob
from ROOT import *
output_filename = glob.glob('*.root')[0]
fin = TFile(output_filename, 'read')
tree = fin.Get('TreeMaker2/PreSelection')
if tree.GetBranch('tracks'):
    print 'File OK'
else:
    print 'no tracks collection, deleting output file'
    os.system('rm ' + output_filename)
    sys.exit(919191)
" > check.py
python check.py
if [[ $? -ne 0 ]]; then
	# this is the exit code for missing branch "tracks" in output file:
    exit 51919
fi

# copy output to eos
echo "gfal-copy output for condor"
if [ -e "/cvmfs/oasis.opensciencegrid.org/mis/osg-wn-client/3.3/current/el6-x86_64/setup.sh" ]; then
    . /cvmfs/oasis.opensciencegrid.org/mis/osg-wn-client/3.3/current/el6-x86_64/setup.sh
fi
for FILE in *.root
do
  echo gfal-copy -n 1 "file:////$PWD/${FILE}" "${OUTDIR}${FILE}"
  gfal-copy -n 1 "file:////$PWD/${FILE}" "${OUTDIR}${FILE}"
  XRDEXIT=$?
  if [[ $XRDEXIT -ne 0 ]]; then
    rm *.root
    echo "exit code $XRDEXIT, failure in xrdcp"
    exit $XRDEXIT
  fi
  rm ${FILE}
done

