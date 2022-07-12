import os,sys

if __name__=="__main__":
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument("path",type=str,default="",help="path containing crab directories")
    parser.add_argument("--grep","-g",type=str,help="grep specific lines")
    args=parser.parse_args()
    
    if args.grep is None:
        lines=os.popen("find "+args.path+" -type d|grep 'crab_projects/[^/]*$'|sort -V").read().split()
    else:
        lines=os.popen("find "+args.path+" -type d|grep 'crab_projects/[^/]*$'|sort -V|grep "+args.grep).read().split()
    if lines[-1]=="": lines=lines[:-1]
    
    for line in lines:
        name=os.path.basename(line)
        progress=os.popen("crab status -d "+line+"|grep -o 'finished.*$'|awk '{print $2}'").read().strip()
        if progress=="": progress="0.0%"
        print name,progress

    
