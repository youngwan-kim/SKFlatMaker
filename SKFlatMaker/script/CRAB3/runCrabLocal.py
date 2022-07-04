import os

def get_unfinished_jobids(args):
    jobids_raw=os.popen("crab status -d "+args.dir+" --long|awk '/ Job State/{flag=1;next;} /^$/{flag=0} {if(flag&&$2!=\"finished\"){print $1}}'").read().split()
    jobids=[]
    for jobid in jobids_raw:
        if jobid.isdigit():
            jobids+=[int(jobid)]
    return jobids

def write_run_script(args):
    script=\
'''#!/bin/bash
export SCRATCH_DIR=${{_CONDOR_SCRATCH_DIR:-$PWD}}
(
    cd {cmssw_base}/src/SKFlatMaker
    source setup.sh
    crab preparelocal -d {crabdir} --destdir $SCRATCH_DIR --proxy $SCRATCH_DIR/x509up_u{uid}
)
cd $SCRATCH_DIR
sh run_job.sh $1
EXITCODE=$?
if [ "$EXITCODE" != 0 ]; then
  echo "Non-zero exit code $EXITCODE" 1>&2
  exit $EXITCODE
fi
ls -1 *.root|sed "p;s/.root/_${{1}}.root/"|xargs -n2 mv
mv *.root {crabdir}/condor
'''.format(cmssw_base=args.cmssw_base,crabdir=args.dir,uid=os.getuid())
    with open(args.dir+"/condor/run.sh","w") as f:
        f.write(script)
    return
    
def write_job_script(args):
    script=\
'''jobbatchname={jobname}
executable=run.sh
output=job$(arguments).out
error=job$(arguments).err
log=job$(arguments).log
request_memory=5000
should_transfer_files=yes
transfer_output_files=""
max_retries=2
on_exit_hold=exitcode!=0
getenv=false
x509userproxy=/tmp/x509up_u{uid}
queue arguments from (
{jobids}
)
'''.format(jobname=os.path.basename(args.dir),jobids="\n".join([str(i) for i in args.jobids]),uid=os.getuid())
    with open(args.dir+"/condor/condor.jds","w") as f:
        f.write(script)
    return

if __name__=="__main__":
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument("--dir","-d",required=True,help="path of crab directory")
    parser.add_argument("--jobids",help="job IDs to submit")
    args=parser.parse_args()

    hostname=os.popen("hostname").read().strip()
    if "tamsa" not in hostname:
        print "This script only runs at tamsa"
        exit(1)

    args.dir=os.path.abspath(args.dir)
    
    if "CMSSW_BASE" not in os.environ:
        print "Please setup cmssw environment"
        exit(1)
    else:
        args.cmssw_base=os.environ["CMSSW_BASE"]

    if args.jobids:
        args.jobids=[int(jobid) for jobid in args.jobids.split(",")]
    else:
        args.jobids=get_unfinished_jobids(args)
        
    if len(args.jobids)==0:
        print "No unfinished job"
        exit(1)

    print "Jobs to run:",args.jobids

    workingdir=args.dir+"/condor"
    if not os.path.exists(workingdir):
        os.makedirs(workingdir)
    write_run_script(args)
    write_job_script(args)
    os.system("cd {};condor_submit condor.jds".format(workingdir))

