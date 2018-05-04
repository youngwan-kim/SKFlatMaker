import os

#txtfilename = 'samplelist_MC.txt'
txtfilename = 'samplelist_DATA_DoubleMuon.txt'
#txtfilename = 'samplelist_DATA_DoubleEG.txt'

SKFlatTag = os.environ['SKFlatTag']

lines = open(txtfilename).readlines()
pdfidshifts = open('PDFInfo.txt').readlines()

isData = "MC"
if ("DATA" in txtfilename) or ("Data" in txtfilename) or ("data" in txtfilename):
  isData = "DATA"
if ("Private" in txtfilename) or ("private" in txtfilename):
  isData = "PrivateMC"

base_dir = SKFlatTag+'/crab_submission_'+isData+'/'

os.system('mkdir -p '+base_dir)
os.system('cp RunSKFlatMaker.py '+base_dir)

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
  pdfshift = "0"
  order = "NLO"
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

  outname = base_dir+'/SubmitCrab__'+sample+'__'+confs+'.py'
  out = open(outname,'w')

  for sk_line in sk_lines:
    if "config.General.requestName" in sk_line:
      if isData=="DATA":
        out.write("config.General.requestName = '"+sample+"__"+confs+"'\n")
      else:
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
    elif 'config.Data.unitsPerJob' in sk_line:
      if isData=="DATA":
        out.write("config.Data.unitsPerJob = 100\n")
      else:
        out.write("config.Data.unitsPerJob = 1\n")
    elif 'config.Data.outputDatasetTag' in sk_line:
      if isData=="DATA":
        period = ""
        if "2017B" in confs:
          period = "B"
        elif "2017C" in confs:
          period = "C"
        elif "2017D" in confs:
          period = "D"
        elif "2017E" in confs:
          period = "E"
        elif "2017F" in confs:
          period = "F"
        out.write("config.Data.outputDatasetTag = 'SKFlat_"+SKFlatTag+"_period"+period+"'\n")
      else:
        out.write("config.Data.outputDatasetTag = 'SKFlat_"+SKFlatTag+"'\n")
    else:
      out.write(sk_line)

  if isData=="DATA":
    out.write("config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/Final/Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt'\n")

  if isData=="PrivateMC":
    out.write("config.Data.inputDBS = 'phys03'\n")

  cmd = 'crab submit -c SubmitCrab__'+sample+'__'+confs+'.py'
  print cmd

  out.close()
