import os,sys

def queryStatus(dasname):
    return os.popen('dasgoclient --query "{} | grep dataset.status"|sort|uniq|tail -n1'.format(dasname)).read().strip()
def queryEvents(dasname):
    return os.popen('dasgoclient --query "{} | grep dataset.nevents"|sort|uniq|tail -n1'.format(dasname)).read().strip()
def queryCompletedAndTotalEvents(dasname):
    words=os.popen('dasgoclient --query "mcm dataset={} | grep mcm.completed_events,mcm.total_events"|sort|uniq|tail -n1'.format(dasname)).read().strip().split()
    if len(words)!=2:
        print("Fail to query {} from McM : {}".format(dasname," ".join(words)))
        return "0/0(0%)"
    return "{}/{}({}%)".format(words[0],words[1],100*int(words[0])/int(words[1]))

if __name__=="__main__":
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument("input",help="DAS name or file with DAS names")
    parser.add_argument("-w","--write",action="store_true",help="overwrite the file instead of printing")
    
    args=parser.parse_args()

    args.inputIsFile=False
    if os.path.exists(args.input):
        args.inputIsFile=True
        with open(args.input) as f:
            lines=f.read().split("\n")
            if lines[-1]=="": lines=lines[:-1]
    elif args.input.count("/")==3:
        lines=[args.input]
    
    out=None
    if args.write and args.inputIsFile:
        out=open(".checkDAS_{}".format(args.input),"w")

    for line_with_comment in lines:
        fields=line_with_comment.split("#",1)
        if len(fields)==2:
            line=fields[0]
            comment="#"+fields[1]
        else:
            line=fields[0]
            comment=""
            
        dasname=line.split()[0] if len(line.split()) else ""
        if dasname.count("/")!=3:
            sys.stdout.write(line+comment+"\n")
            if out: out.write(line+comment+"\n")
        else:
            status=queryStatus(dasname)
            if status=="VALID":
                events=queryEvents(dasname)
            elif status=="PRODUCTION":
                events=queryCompletedAndTotalEvents(dasname)
            else:
                print("Unknown status {} for {}".format(status,dasname))
                events=status
            sys.stdout.write("{}\t{}{}\n".format(dasname,events,comment))
            if out: out.write("{}\t{}{}\n".format(dasname,events,comment))
    if out:
        out.close()
        os.system("mv {} {}".format(".checkDAS_{}".format(args.input),args.input))
    exit(0)
