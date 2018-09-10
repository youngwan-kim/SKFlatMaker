import os

SKFlatTag = os.environ['SKFlatTag']

os.system('ls *.log &> tmp.txt')

lines = open('tmp.txt').readlines()

out=open('SampleLogs.html','w')
out.write('<table border = 1>\n')
out.write('  <tr>\n')
out.write('    <th><a name="LHE Log files">LHE Log files</a></th>\n')
out.write('  </tr>\n')
out.write('  <tr>\n')
out.write('    <th><a href="https://docs.google.com/spreadsheets/d/17tNM-1pYh8ktM5cH_pFI8oVXuYR6QzPCtY4YWzLdMKY/edit#gid=0">Spreadsheet Link</a></th>\n')
out.write('  </tr>\n')

for line in lines:
  line = line.strip('\n')
  samplenameonly = line.replace('.log','')

  out.write('<tr>\n')
  out.write('  <td><a href="'+line+'">'+samplenameonly+'</td>\n')
  out.write('</tr>\n')

  urllink = 'http://jskim.web.cern.ch/jskim/SKFlat/LHELogs/'+SKFlatTag+'/'+line
  print '=HYPERLINK("'+urllink+'","'+samplenameonly+'")'

out.close()

print 'http://jskim.web.cern.ch/jskim/SKFlat/LHELogs/'+SKFlatTag+'/SampleLogs.html'

os.system('mkdir -p /eos/user/j/jskim/www/SKFlat/LHELogs/'+SKFlatTag)
os.system('cp * /eos/user/j/jskim/www/SKFlat/LHELogs/'+SKFlatTag)
