import os

lines = open('samplelist.txt').readlines()
pdfidshifts = open('PDFIDShift.txt').readlines()

hostname = os.environ['HOSTNAME']

os.system('mkdir -p logs')

for line in lines:

  if "#" in line:
    continue

  line = line.strip('\n')
  samplePDs = line.split("/")

  sample = samplePDs[1]

  dascmd = '/cvmfs/cms.cern.ch/common/dasgoclient --limit=1 --query="file dataset='+line+'"'
  print dascmd

  os.system(dascmd+' > tmp.txt')

  filepaths = open('tmp.txt').readlines()

  if len(filepaths)==0:
    print "No file for : "+sample
    break

  filepath = filepaths[0].replace('/xrootd','').strip('\n')

  if 'lxplus' not in hostname:
    filepath = 'root://cms-xrd-global.cern.ch/'+filepath

  cmd = 'cmsRun RunSKFlatMaker.py '+filepath

  for shiftsamples in pdfidshifts:
    words = shiftsamples.split()
    if words[0]==sample:
      this_shift = words[1]
      if "?" in this_shift:
        this_shift = "0"
      cmd = cmd+' '+this_shift
      break

  final_cmd = cmd+' > logs/'+sample+'.log'
  print final_cmd
  os.system(final_cmd)

os.system('rm tmp.txt')
