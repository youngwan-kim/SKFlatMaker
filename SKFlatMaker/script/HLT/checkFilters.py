import os,sys,argparse
parser=argparse.ArgumentParser()
parser.add_argument("path",type=str)
parser.add_argument("--era",default="2016,2017,2018")
parser.add_argument("--debug",action="store_true")
args=parser.parse_args()

import hlt2016,hlt2017,hlt2018
paths=[]
if "2016" in args.era:
    paths+=[("2016",x,getattr(hlt2016.fragment,x)) for x in dir(hlt2016.fragment) if args.path in x]
if "2017" in args.era:
    paths+=[("2017",x,getattr(hlt2017.fragment,x)) for x in dir(hlt2017.fragment) if args.path in x]
if "2018" in args.era:
    paths+=[("2018",x,getattr(hlt2018.fragment,x)) for x in dir(hlt2018.fragment) if args.path in x]
for p in paths:
    if args.debug: print p
    print p[0], p[1]
    print [x for x in str(p[2]).split("+") if "Filter" in x]
