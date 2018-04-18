import os

txtfilename = 'samplelist_MC.txt'
lines = open(txtfilename).readlines()
pdfidshifts = open('PDFInfo.txt').readlines()

isData = "MC"
if ("DATA" in txtfilename) or ("Data" in txtfilename) or ("data" in txtfilename):
  isData = "DATA"
if ("Private" in txtfilename) or ("private" in txtfilename):
  isData = "PrivateMC"

os.system('cp RunSKFlatMaker.py crab_submission_'+isData+'/')

hostname = os.environ['HOSTNAME']

for line in lines:

  if "#" in line:
    continue

  line = line.strip('\n')
  samplePDs = line.split("/")

  sample = samplePDs[1]
  confs = samplePDs[2]

  # cmsRun RunSKFlatMaker.py <DATA/MC/PrivateMC> <PDFIDShift;M1> <NLO/LO> <powheg/madgraph0/madgraph1000>

  # DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8 0 LO  madgraph0
  # DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8  0 NLO madgraph1000
  pdfshift = ""
  order = ""
  gentype = ""

  if "MC" in isData:
    for shiftsamples in pdfidshifts:
      words = shiftsamples.split()
      if words[0]==sample:
        pdfshift = words[1]
        order = words[2]
        gentype = words[3]

  if order=="?":
    order="NLO"

  sk_lines = open('skeleton/SubmitCrab.py').readlines()

  #os.system('mkdir -p crab_submission_'+isData+'/'+sample+'/')
  #outname = 'crab_submission_'+isData+'/'+sample+'/'+confs+'.py'

  os.system('mkdir -p crab_submission_'+isData+'/')
  outname = 'crab_submission_'+isData+'/SubmitCrab__'+sample+'__'+confs+'.py'
  out = open(outname,'w')

  for sk_line in sk_lines:
    if "config.General.requestName" in sk_line:
      out.write("config.General.requestName = '"+sample+"'\n")
    elif "config.JobType.pyCfgParams" in sk_line:
      out.write("config.JobType.pyCfgParams = ['sampletype="+isData+"','PDFIDShift="+pdfshift+"','PDFOrder="+order+"','PDFType="+gentype+"']\n")
    elif "config.Data.inputDataset" in sk_line:
      out.write("config.Data.inputDataset = '"+line+"'\n")
    elif 'config.Data.splitting' in sk_line:
      if isData=="DATA":
        out.write("config.Data.splitting = 'LumiBased'\n")
      else:
        out.write("config.Data.splitting = 'FileBased'\n")
    else:
      out.write(sk_line)

  if isData=="PrivateMC":
    out.write("config.Data.inputDBS = 'phys03'\n")

  out.close()
