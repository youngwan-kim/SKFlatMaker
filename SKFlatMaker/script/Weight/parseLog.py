import os,sys,argparse
import re

parser = argparse.ArgumentParser(description="parsing SKFlatMaker log to extract weight INDEXs")
parser.add_argument("log",help="log file or directory")
parser.add_argument("--verbose","-v",action="count",default=0,help="-v verbose mode, -vv debug mode")
args=parser.parse_args()

PrintLevel=args.verbose

class Weight:
    def __init__(self,words=[],isRoot=False):
        self.name=""
        self.PDF=""
        self.memberID=-1
        self.index=-1
        self.muF=-1
        self.muR=-1
        self.alpsfact=-1
        self.combine=""
        self.description=""
        self.subs=[]
        self.fullline=" ".join(words)
        self.PrintLevel=0
        if "PrintLevel" in globals():
            self.PrintLevel=PrintLevel
        if type(words) is str:
            words=words.strip().split()
        self.Parse(words,isRoot)
        
    def __getitem__(self,key):
        return self.subs[key]

    def __len__(self):
        return len(self.subs)

    def Find(self,size=-1,PDF="",muF=-1,muR=-1,name=""):
        out=[]
        if (size==-1 or len(self)==size) and (PDF=="" or self.PDF==str(PDF)) and (muF==-1 or self.muF==muF) and (muR==-1 or self.muR==muR) and (name=="" or self.name==name):
            out+=[self]
        for sub in self.subs:
            out+=sub.Find(size=size,PDF=PDF,muF=muF,muR=muR,name=name)
        return out

    def HasSubgroup(self):
        for sub in self.subs:
            if sub.IsGroup(): return True
        return False

    def IsGroup(self):
        if len(self): return True
        else: return False

    def Parse(self,words=[],isRoot=False):
        if isRoot: level=0
        else: level=-1
        buf=""
        for word in words:
            if word=="[SKFlatMaker::endRun,lines]": continue
            buf+=" "+word
            if "<weight" in word:
                level+=1
                if level==0:
                    buf=""
                    continue
            elif word=="</weight>" or word=="</weightgroup>":
                level-=1
                if level==-1:
                    break;
                elif level==0:
                    self.subs+=[Weight(buf)]
                    buf=""
            if level==0:
                short="''".join('""'.join(buf.split('"')[::2]).split("'")[::2])
                if len(short.split("="))==2 and short.count("'")%2==0 and short.count('"')%2==0 and not isRoot:
                    self.Set(buf)
                    buf=""
            elif level<0:
                buf=""
            
    def Print(self,depth=0,verbose=False):
        line=[]
        if self.name!="": line+=["name:"+self.name]
        if self.index!=-1: line+=["index:"+str(self.index)]
        if self.PDF!="": line+=["PDF:"+self.PDF]
        if self.memberID!=-1: line+=["memberID:"+str(self.memberID)]
        if self.muF!=-1: line+=["muF:"+str(self.muF)]
        if self.muR!=-1: line+=["muR:"+str(self.muR)]
        if self.alpsfact!=-1: line+=["alpsfact:"+str(self.alpsfact)]
        if self.combine!="": line+=["combine:"+self.combine]
        if len(self): line+=["nsub:"+str(len(self))]
        if self.description!="": line+=[self.description]
        print(" "*depth+" ".join(line))
        if self.HasSubgroup() or len(self)<10 or verbose:
            for sub in self.subs:
                sub.Print(depth+2)
        else:
            self[0].Print(depth+2)
            self[1].Print(depth+2)
            print(" "*(depth+2)+"...")
            self[-1].Print(depth+2)
            
    def Set(self,line):
        if "=" in line:
            temp=line.split("=")
            words=[]
            for word in temp:
                if len(words)==0: 
                    words=[word]
                    continue
                if words[-1].count("'")%2==1 or words[-1].count('"')%2==1:
                    words[-1]+="="+word
                else:
                    words+=[word]
            if len(words)!=2:
                print("[Weight::Set] ERROR Wrong line "+line)
                exit(1)
            left=words[0].strip().lower()
            right=words[1].strip().rstrip(">").strip('"').strip("'")
        if left == "muf" or left == "facscfact": self.muF=float(right.replace("d","e"))
        elif left == "mur" or left == "renscfact": self.muR=float(right.replace("d","e"))
        elif left == "alpsfact": self.alpsfact=float(right.replace("d","e"))
        elif left == "combine": self.combine=right
        elif left == "name": self.name=right
        elif left == "id": self.index=int(right)
        elif left == "memberid": self.memberID=int(right)
        elif left == "pdf" or left == "lhapdf": 
            if self.PDF!="": self.description+=" "+line                
            else: self.PDF=right
        else:
            self.description+=" "+line
            if self.PrintLevel>1:
                print("[Weight::Set] WARNING unknown key "+words[0].strip())
            


logs=[]
if os.path.isfile(args.log):
    logs=[args.log]
elif os.path.isdir(args.log):
    logs=[os.path.join(args.log,f) for f in os.listdir(args.log) if os.path.isfile(os.path.join(args.log,f))]
else:
    raise OSError(2, 'No such file or directory', args.log)

for log in logs:
    with open(log) as f:
        lines=f.readlines()
        words=[]
        for line in lines:
            words+=line.split()

    w=Weight(words,True)
    if PrintLevel>0:
        w.Print()
    
    scale=[]
    alpsfact=[]
    warning=""
    combine=""
    pdf=[]
    alphaS=[]
    PSSyst=[]
    asScale=0
    lhaid=""

    ###############
    #### LHEID ####
    for line in lines:
        m=re.search(r"([0-9]+)[ \t]*=[ \t]*lhaid",line)
        if m:
            lhaid=m.group(1)
            break
        m=re.search(r"lhans1[ \t]+([0-9]+)",line)
        if m:
            lhaid=m.group(1)
            break


    ###############
    #### asScale ##
    if lhaid!="":
        lhapdf=os.popen("grep "+lhaid+" $LHAPDF_DATA_PATH/pdfsets.index|awk '{print $2}'").read().strip()
        lhapdfinfo=os.popen("head -n1 $LHAPDF_DATA_PATH/"+lhapdf+"/"+lhapdf+".info").read().strip()
        if "hessian" in lhapdfinfo.lower() or "hessian" in lhapdf.lower(): combine="hessian"
        else: combine="gaussian"
        if "0.116" in lhapdfinfo and "0.120" in lhapdfinfo:
            asScale=0.75
        elif "0.117" in lhapdfinfo and "0.119" in lhapdfinfo:
            asScale=1.5
        if PrintLevel>0:
            print("---- Default PDF, "+lhapdf+" "+lhaid+" -----------")
            print lhapdfinfo
    else:
        print("> Cannot find default lhaid.. assume lhaid=306000")
        warning+=" Assume lhaid=306000."
        lhaid="306000"

    
    ###############
    #### Scale ####
    if PrintLevel>0:
        print("> Search candidates for scale variation group")
    cands=w.Find(size=9)
    if len(cands)>1:
        if PrintLevel>0:
            print("> Multiple candidates for scale variation group")
            for cand in cands:
                cand.Print()
            print("------------------------------")
        newcands=[]
        for cand in cands:
            if len(cand.Find(PDF=int(lhaid))):
                newcands+=[cand]
        cands=newcands
        if PrintLevel>0:
            print("> Selected scale variation group")
            for cand in cands:
                cand.Print()
            print("------------------------------")

    if len(cands)==0:
        if PrintLevel>0:
            print("> No candidates for scale variation group.. use method 2")
        cand=Weight()
        cand.subs+=w.Find(muF=1.0,muR=1.0,PDF=lhaid)[:1]
        cand.subs+=w.Find(muF=1.0,muR=2.0,PDF=lhaid)[:1]
        cand.subs+=w.Find(muF=1.0,muR=0.5,PDF=lhaid)[:1]
        cand.subs+=w.Find(muF=2.0,muR=1.0,PDF=lhaid)[:1]
        cand.subs+=w.Find(muF=2.0,muR=2.0,PDF=lhaid)[:1]
        cand.subs+=w.Find(muF=2.0,muR=0.5,PDF=lhaid)[:1]
        cand.subs+=w.Find(muF=0.5,muR=1.0,PDF=lhaid)[:1]
        cand.subs+=w.Find(muF=0.5,muR=2.0,PDF=lhaid)[:1]
        cand.subs+=w.Find(muF=0.5,muR=0.5,PDF=lhaid)[:1]
        #cand.Print(verbose=True)
        if len(cand)==9:
            cands+=[cand]
        
    if len(cands)==1:
        scale=[cands[0][i].index for i in range(9)]

    ###############
    #### Scale ####
    if PrintLevel>0:
        print("> Search candidates for alpsfact variation group")
    cands=w.Find(size=2)
    newcands=[]
    for cand in cands:
        if cand.subs[0].alpsfact!=cand.subs[1].alpsfact:
            newcands+=[cand]
    cands=newcands
    if len(cands)==1:
        alpsfact=[cands[0][i].index for i in range(2)]

    ###############
    #### PDF ######
    if PrintLevel>0:
        print("> Search candidates for PDF error set group")
    cands=w.Find(size=100)
    cands+=w.Find(size=101)
    cands+=w.Find(size=102)
    cands+=w.Find(size=103)    
    if len(cands)>1:
        newcands=[]
        for cand in cands:
            if len(cand.Find(PDF=int(lhaid)+1)) and len(cand.Find(PDF=int(lhaid)+100)):
                newcands+=[cand]
        cands=newcands

    if len(cands)==0:
        cand=Weight()
        for i in range(1,103):
            cand.subs+=w.Find(PDF=int(lhaid)+i)
        if len(cand)>99:
            cands+=[cand]

    if len(cands)==0:
        altlhaids=[]
        if lhaid=="306000": altlhaids=["305800"]
        for altlhaid in altlhaids:
            print "> Cannot find pdf variations for {}. Try {}.".format(lhaid,altlhaid)
            cand=Weight()
            for i in range(101):
                cand.subs+=w.Find(PDF=int(altlhaid)+i)
            if len(cand)>99:
                cands+=[cand]
                warning+=" Use {} for the PDF variation set.".format(altlhaid)
                break
            
    if len(cands)==1: 
        final=cands[0]
        #final.Print()
        if len(final)==100:
            pdf=[final[i].index for i in range(0,100)]
        elif len(final)==101:
            pdf=[final[i].index for i in range(1,101)]
        elif len(final)==102:
            pdf=[final[i].index for i in range(0,100)]
            alphaS=[final[-2].index,final[-1].index]
        elif len(final)==103:
            pdf=[final[i].index for i in range(1,101)]
            alphaS=[final[-2].index,final[-1].index]
    
    if len(alphaS)==0:
        if lhaid=="306000":
            alphaSA=w.Find(PDF="323300")
            alphaSB=w.Find(PDF="323500")
            if len(alphaSA) and len(alphaSB):
                alphaS=[alphaSA[0].index,alphaSB[0].index]
                asScale=1.5

    ##### Pythia parton shower systematic ####
    cands=w.Find(name="fsr:murfac=0.5")+w.Find(name="fsr:murfac=2.0")+w.Find(name="isr:murfac=0.5")+w.Find(name="isr:murfac=2.0")
    PSSyst=[cand.index for cand in cands]
    if PSSyst!=[4,5,26,27]:
        print "> Unusual Pythia weights("+PSSyst.__str__()+"). Please Check."

    ###############
    #### output ###
    out=["## lhaid="+lhaid+" "+combine+warning]
    if len(scale)==9:
        if scale[-1]-scale[0]==9-1:
            out+=["Scale\t{},{}".format(scale[0],scale[-1])]
        else:
            out+=["Scale\t"+",".join([str(i) for i in scale])]

    if len(alpsfact):
        if alpsfact[-1]-alpsfact[0]==len(alpsfact)-1:
            out+=["alpsfact\t{},{}".format(alpsfact[0],alpsfact[-1])]
        else:
            out+=["alpsfact\t"+",".join([str(i) for i in alpsfact])]

    if len(pdf)>1: 
        if pdf[-1]-pdf[0]==len(pdf)-1:
            out+=["PDF\t{},{}".format(pdf[0],pdf[-1])]
        else:
            out+=["PDF\t"+",".join([str(i) for i in pdf])]

    if len(alphaS)>1: 
        if alphaS[-1]-alphaS[0]==len(alphaS)-1:
            out+=["AlphaS\t{},{}".format(alphaS[0],alphaS[-1])]
        else:
            out+=["AlphaS\t"+",".join([str(i) for i in alphaS])]

    if asScale!=0:
        out+=["AlphaSScale\t{}".format(asScale)]

    if len(PSSyst)>1:
        if PSSyst[-1]-PSSyst[0]==len(PSSyst)-1:
            out+=["PSSyst\t{},{}".format(PSSyst[0],PSSyst[-1])]
        else:
            out+=["PSSyst\t"+",".join([str(i) for i in PSSyst])]

    ## MiNNLO specific
    cands=w.Find(name="sthw2_variation")
    if len(cands):
        out+=["sthw2\t"+",".join([str(c.index) for c in cands[0].subs])]
    cands=w.Find(name="q0_variation")
    if len(cands):
        out+=["q0\t"+",".join([str(c.index) for c in cands[0].subs])]
    cands=w.Find(name="largeptscales_variation")
    if len(cands):
        out+=["largeptscales\t"+",".join([str(c.index) for c in cands[0].subs])]

    outfile=os.environ["SKFlatWD"]+"SKFlatMaker/script/Weight/data/"+os.path.basename(log.replace(".log",".txt"))
    if os.path.exists(outfile):
        print "-- The file already exists. "+outfile
        print "--- diff old new"
        diff=os.system("bash -c \"diff "+outfile+" <( "+";".join(["echo '"+l+"'" for l in out])+" )\"")
        if diff==0:
            print "--- No change needed. Continue..."
            continue
    print "=== Output for {} ===".format("/"+os.path.basename(log).replace(".log","__MINIAODSIM").replace("__","/"))
    for line in out: print line
    yes=raw_input("---> "+outfile+"\n(y/n):")
    if yes.lower()=="y":
        os.system("mkdir -p "+os.path.dirname(outfile))
        with open(outfile,"w") as f:
            for line in out:
                f.write(line+"\n")
    
