import os, json

def get_unfinished_jobids(args):
    jobids_raw=os.popen("crab status -d "+args.dir+" --long|awk '/ Job State/{flag=1;next;} /^$/{flag=0} {if(flag&&$2!=\"finished\"){print $1}}'").read().split()
    jobids=[]
    for jobid in jobids_raw:
        if jobid.isdigit():
            jobids+=[int(jobid)]
    return jobids

def fetch_file(filename):
    das=json.loads(os.popen("dasgoclient --query 'site file={}|grep site.rses'|tail -n1".format(filename)).read())
    target="/gv0/Users/{user}/runCrabLocal{path}".format(user=os.getlogin(),path=filename)
    adler32_remote=os.popen("dasgoclient --query 'site file={} |grep site.adler32' | sort | tail -n1".format(filename)).read().strip()
    print "> Fetching",filename
    if os.path.exists(target):
        adler32_local=os.popen("xrdadler32 "+target+" |awk '{print $1}'").read().strip()
        if adler32_local==adler32_remote:
            print "File exists with correct checksum"
            return 0
        else:
            print "File with wrong checksum {}:{}. Remove {}".format(adler32_remote,adler32_local,target)
            os.system("rm "+target)
    for site in das.keys():
        source=das[site][0]
        #if "ucsd" in source: continue ## probably corrupted
        if not os.path.exists(os.path.dirname(target)):
            os.makedirs(os.path.dirname(target))
        env="PYTHONHOME=/usr/lib64/python2.7 PYTHONPATH=/usr/lib64/python2.7/site-packages:/usr/lib64/python2.7:/usr/lib/python2.7/site-packages:/usr/lib64/python2.7/lib-dynload LD_LIBRARY_PATH=/usr/lib64/gfal2-plugins"
        cmd="gfal-copy {source} file://{target}".format(source=source,target=target)
        print cmd
        rt=os.system(env+" "+cmd)
        if rt!=0:
            print "Fail to run",cmd
            if os.path.exists(target):
                os.system("rm "+target)
        if rt==0:
            print "Calculating checksum..."
            adler32_local=os.popen("xrdadler32 "+target+" |awk '{print $1}'").read().strip()
            if adler32_local==adler32_remote:
                return rt
            else:
                print "File with wrong checksum {}:{}. Remove {}".format(adler32_remote,adler32_local,target)
                print "Try next source"
                rt=-1
    return rt

if __name__=="__main__":
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument("--dir","-d",help="path of crab directory")
    parser.add_argument("--file","-f",help="path of file to fetch")
    parser.add_argument("--jobids",help="job IDs to submit")
    args=parser.parse_args()

    hostname=os.popen("hostname").read().strip()
    if "tamsa" not in hostname:
        print "This script only runs at tamsa"
        exit(1)

    if "CMSSW_BASE" not in os.environ:
        print "Please setup cmssw environment"
        exit(1)
    else:
        args.cmssw_base=os.environ["CMSSW_BASE"]

    if args.dir:
        args.dir=os.path.abspath(args.dir)
        if args.jobids:
            args.jobids=[int(jobid) for jobid in args.jobids.split(",")]
        else:
            args.jobids=get_unfinished_jobids(args)
        
        if len(args.jobids)==0:
            print "No unfinished job"
            exit(1)

        os.system("crab preparelocal -d {crabdir} --proxy /tmp/x509up_u{uid} &>/dev/null ;cd {crabdir}/local; tar zxf input_files.tar.gz &>/dev/null".format(crabdir=args.dir,uid=os.getuid()))

        for jobid in args.jobids:
            with open("{}/local/job_input_file_list_{}.txt".format(args.dir,jobid)) as f:
                line=f.read()
                for char in ['"',"[","]",","]:
                    line=line.replace(char," ")
                files=line.split()
            print ">> Job",jobid
            for filename in files:
                fetch_file(filename)
    elif args.file:
        fetch_file(args.file)

    else:
        print "Neither dir or file argument provided"
        exit(1)
