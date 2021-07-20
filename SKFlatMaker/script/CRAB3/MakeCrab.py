import os, sys

txtfilename = sys.argv[1]
#txtfilename = '2016_DATA.txt'
#txtfilename = '2016_MC.txt'
#txtfilename = '2017_DATA.txt'
#txtfilename = '2017_MC.txt'
#txtfilename = '2018_DATA.txt'
#txtfilename = '2018_MC.txt'

SKFlatTag = os.environ['SKFlatTag']
SKFlatWD = os.environ['SKFlatWD']

era = "-1"
if "2016preVFP" in txtfilename:
  era = "2016preVFP"
elif "2016postVFP" in txtfilename:
  era = "2016postVFP"
elif "2017" in txtfilename:
  era = "2017"
elif "2018" in txtfilename:
  era = "2018"
else:
  print "Wrong era from txtfilename : "+txtfilename
  sys.exit()
year=era[0:4]

lines = open(txtfilename).readlines()

str_sample = "DATA"
isData = False
isPrivateMC = False
if ("DATA" in txtfilename) or ("Data" in txtfilename) or ("data" in txtfilename):
  str_sample = "DATA"
  isData = True
  isPrivateMC = False
else:
  str_sample = "MC"
  isData = False
  isPrivateMC = False
if ("Private" in txtfilename) or ("private" in txtfilename):
  str_sample = "PrivateMC"
  isData = False
  isPrivateMC = True

base_dir = SKFlatTag+'/'+era+'/crab_submission_'+str_sample+'/'

os.system('mkdir -p '+base_dir)
os.system('cp '+SKFlatWD+'/SKFlatMaker/test/RunSKFlatMaker.py '+base_dir)

hostname = os.environ['HOSTNAME']

for line in lines:

  if "#" in line:
    continue
  if line.count('/')!=3:
    continue

  line = line.strip('\n')

  #### incase the line contains # of event,
  line = line.split()[0]

  samplePDs = line.split("/")

  #continue if blank line or invalid format
  if len(samplePDs)<4:
    continue
  
  sample = samplePDs[1]
  confs = samplePDs[2]

  #### get extension info
  #### _ext1-v2
  extension = ''
  for w in range(0,len(confs)):
    if confs[w:w+4]=='_ext':
      extension = confs[w:]
      break

  cmd = 'crab submit -c SubmitCrab__'+sample+'__'+confs+'.py'

  ArgsListString = "['sampletype="+str_sample+"','era="+era+"'"
  # ['sampletype="+str_sample+"','PDFIDShift="+pdfshift+"','PDFOrder="+order+"','PDFType="+gentype+"','era="+era+"']
  if not isData:
    HasFile = os.path.isfile(SKFlatWD+'/SKFlatMaker/script/MCPDFInfo/'+era+'/'+sample+'.txt')

    if HasFile:

      MCInfoLines = open(SKFlatWD+'/SKFlatMaker/script/MCPDFInfo/'+era+'/'+sample+'.txt').readlines()
      words = MCInfoLines[0].split()
      # 1001,1009 gaussian  2001,2100 2101,2102 1.5,1.5
      ScaleIDRange = words[0]
      PDFErrorType = words[1]
      PDFErrorIDRange = words[2]
      PDFAlphaSIDRange = words[3]
      PDFAlphaSScaleValue = words[4]

      #ArgsListString += ",'PDFErrorType="+PDFErrorType+"'"
      ArgsListString += ",'ScaleIDRange="+ScaleIDRange+"'"
      ArgsListString += ",'PDFErrorIDRange="+PDFErrorIDRange+"'"
      ArgsListString += ",'PDFAlphaSIDRange="+PDFAlphaSIDRange+"'"
      ArgsListString += ",'PDFAlphaSScaleValue="+PDFAlphaSScaleValue+"'"

      AdditionalWeightList=[]
      for i in range(1,len(MCInfoLines)):
        AdditionalWeightList += ["{}[{}]".format(MCInfoLines[i].split()[0],MCInfoLines[i].split()[1])]
      if len(AdditionalWeightList):
        ArgsListString += ",'AdditionalWeights="+",".join(AdditionalWeightList)+"'"

      if "999" in MCInfoLines[0]:
        print '#### Has Issue : '+sample
        #continue

    #### If WR sample
    elif 'WRtoNLtoLLJJ' in sample:
      ArgsListString += ",'ScaleIDRange=1001,1045'"
      ArgsListString += ",'PDFErrorIDRange=1046,1146'"
      ArgsListString += ",'PDFAlphaSIDRange=1147,1148'"
      ArgsListString += ",'PDFAlphaSScaleValue=0.75,0.75'"

    else:
      ### to avoid exit()
      print '#### No MCPDFInfo : '+sample
      #ArgsListString += ",'PDFErrorType="+"hessian"+"'"

  ArgsListString += "]"


  sk_lines = open('skeleton/SubmitCrab.py').readlines()

  outname = base_dir+'/SubmitCrab__'+sample+'__'+confs+'.py'
  out = open(outname,'w')

  for sk_line in sk_lines:
    if "config.General.requestName" in sk_line:
      if isData:
        out.write("config.General.requestName = '"+sample+"__"+confs+"'\n")
      else:
        out.write("config.General.requestName = '"+sample+extension+"'\n")
    elif "config.JobType.pyCfgParams" in sk_line:
      out.write("config.JobType.pyCfgParams = "+ArgsListString+"\n")
    elif "config.Data.inputDataset" in sk_line:
      out.write("config.Data.inputDataset = '"+line+"'\n")
      if isPrivateMC:
        out.write("config.Data.inputDBS = 'phys03'\n")
        out.write("config.Data.ignoreLocality = True\n")
        out.write("config.Site.whitelist = ['T2*']\n")

    elif 'config.Data.splitting' in sk_line:
      if isData:
        out.write("config.Data.splitting = 'LumiBased'\n")
      else:
        out.write("config.Data.splitting = 'FileBased'\n")
    elif 'config.Data.unitsPerJob' in sk_line:
      if isData:
        out.write("config.Data.unitsPerJob = 50\n")
      else:
        out.write("config.Data.unitsPerJob = 1\n")
    elif 'config.Data.outputDatasetTag' in sk_line:
      if isData:
        period = ""
        if year+"A" in confs:
          period = "A"
        elif year+"B" in confs:
          period = "B"
          if year=="2016":
            if "ver1" in confs:
              period = "B_ver1"
            if "ver2" in confs:
              period = "B_ver2"
        elif year+"C" in confs:
          period = "C"
        elif year+"D" in confs:
          period = "D"
        elif year+"E" in confs:
          period = "E"
        elif year+"F" in confs:
          period = "F"
        elif year+"G" in confs:
          period = "G"
        elif year+"H" in confs:
          period = "H"

        out.write("config.Data.outputDatasetTag = 'SKFlat_"+SKFlatTag+"_period"+period+"'\n")
      else:
        out.write("config.Data.outputDatasetTag = 'SKFlat_"+SKFlatTag+"'\n")

    elif 'config.Data.outLFNDirBase' in sk_line:
      out.write("config.Data.outLFNDirBase = '/store/user/%s/SKFlat/"+era+"/' % (getUsername())\n")

    elif 'config.Site.storageSite' in sk_line:
      out.write("config.Site.storageSite = 'T3_KR_KNU'\n")
    else:
      out.write(sk_line)

  if isData:
    if era=="2016preVFP":
      out.write("config.Data.lumiMask = 'https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions16/13TeV/Legacy_2016/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt'\n")
    elif era=="2016postVFP":
      out.write("config.Data.lumiMask = 'https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions16/13TeV/Legacy_2016/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt'\n")
    elif era=="2017":
      #### https://twiki.cern.ch/twiki/bin/view/CMS/PdmVLegacy2017Analysis
      out.write("config.Data.lumiMask = 'https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions17/13TeV/Legacy_2017/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt'\n")
    elif era=="2018":
      #### https://twiki.cern.ch/twiki/bin/view/CMS/PdmVLegacy2018Analysis
      out.write("config.Data.lumiMask = 'https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions18/13TeV/Legacy_2018/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt'\n")
    else:
      print "Wrong era : "+era

  print cmd

  out.close()
