import os,sys,argparse

def getLog(dasname,era=None,isPrivate=False,overwrite=False,fileindex=0):
  if dasname.count('/')!=3:
    print "Wrong DAS name format"
    return False

  dasname_ = "__".join(dasname.split("/")[1:3])
  outfilepath="logs/"+dasname_+".log"
  if os.path.exists(outfilepath):
    if os.system("tail -n1 "+outfilepath+"|grep -q SKFlatMaker::endJob")==0 and not overwrite:
      print "Valid log file exists for {}. Skip...".format(dasname)
      return True

  if era==None:
    if "UL16" in dasname and "APV" in dasname:
      era = "2016preVFP"
    elif "UL16" in dasname:
      era = "2016postVFP"
    elif "UL17" in dasname:
      era = "2017"
    elif "UL18" in dasname:
      era = "2018"
  if era==None:
    dasout=json.load(os.popen('dasgoclient --query="mcm dataset='+dasname+'|grep mcm.sequences"'))
    mcmera=dasout[0]["era"]
    if mcmera=="Run2_2016_HIPM":
      era="2016preVFP"
    elif mcmera=="Run2_2016":
      era="2016postVFP"
    elif mcmera=="Run2_2017":
      era="2017"
    elif mcmera=="Run2_2018":
      era="2018"
  if era==None:
    print "Cannot determine era. Please use a filename with specific era."
    return False

  SKFlatWD = os.environ['SKFlatWD']

  dascmd = 'dasgoclient --limit=10 --query="file dataset='+dasname+'"'
  if isPrivate:
    dascmd = 'dasgoclient --limit=10 --query="file dataset='+dasname+' instance=prod/phys03"'
  print dascmd

  filepaths=os.popen(dascmd).readlines()
  if len(filepaths)==0:
    print "No file for : "+dasname
    return False

  filepath = filepaths[fileindex].replace('/xrootd','').strip('\n')
  filepath = 'root://cms-xrd-global.cern.ch/'+filepath

  cmd = 'cmsRun '+SKFlatWD+'/SKFlatMaker/test/RunSKFlatMaker.py maxEvents=1 inputFiles='+filepath

  weightmap = SKFlatWD+'/SKFlatMaker/script/Weight/data/'+dasname_+'.txt'
  if os.path.isfile(weightmap):
    cmd = cmd+' weightmap='+weightmap

  if isPrivate:
    cmd = cmd+' sampletype=PrivateMC'
  else:
    cmd = cmd+' sampletype=MC'

  cmd = cmd+' era='+era
  
  final_cmd = cmd+' > '+outfilepath
  print final_cmd
  os.system('mkdir -p logs')
  ret=os.system(final_cmd)
  if ret!=0:
    return -1
  return True

if __name__=="__main__":
  parser=argparse.ArgumentParser()
  parser.add_argument("var",nargs="+",help="DAS name or file name of sample list")
  parser.add_argument("--force","-f",action="store_true",help="overwrite existing logs")
  args=parser.parse_args()

  for arg in args.var:
    if os.path.exists(arg):
      lines = open(arg).readlines()
      era = None
      if "2016preVFP" in arg or "2016a" in arg: 
        era = "2016preVFP"
      elif "2016postVFP" in arg or "2016b" in arg: 
        era = "2016postVFP"
      elif "2017" in arg: 
        era = "2017"
      elif "2018" in arg: 
        era = "2018"
      
      for line in lines:
        line=line.split("#")[0]
        line=line.strip()
        if line=="": continue
        line=line.split()[0]
        if line.count('/')!=3:
          continue
        for i in range(3):
          if getLog(line,era=era,isPrivate=("Private" in arg),overwrite=args.force,fileindex=i)!=-1:
            break;

    elif arg.count('/')==3:
      for i in range(3):
        if getLog(arg,overwrite=args.force,fileindex=i)!=-1:
          break;

    else:
      print "Wrong argument "+arg+". No such file or wrong DAS name format"
      continue
    
