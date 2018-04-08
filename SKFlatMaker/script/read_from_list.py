import os

#lines = open('varlist_Muon.txt').readlines()
#lines = open('varlist_Electron.txt').readlines()
#lines = open('varlist_Jet.txt').readlines()
#lines = open('varlist_Trigger.txt').readlines()
lines = open('varlist_LHE.txt').readlines()

# 0 = header
# 1 = src, initilize
# 2 = src, Branch
for WhatToPrint in range(0,3):

  for line in lines:
    words = line.split()

    vartype = words[0]
    varname = words[1]
    isvector = words[2]

    if isvector=="1":
      if WhatToPrint==0:
        print '  vector<'+vartype+'> '+varname+';'
      elif WhatToPrint==1:
        print '  '+varname+'.clear();'
      elif WhatToPrint==2:
        print '    DYTree->Branch("'+varname+'", "vector<'+vartype+'>", &'+varname+');'

    else:
      if WhatToPrint==0:
        print '  '+vartype+' '+varname+';'
      elif WhatToPrint==1:
        print '  '+varname+'=-999;'
      elif WhatToPrint==2:

        vartypechar = ""
        if vartype=="double":
          vartypechar = "D"
        elif vartype=="int":
          vartypechar = "I"
        elif vartype=="bool":
          vartypechar = "O"

        print '    DYTree->Branch("'+varname+'", &'+varname+', "'+varname+'/'+vartypechar+'");'
  print '\n\n\n'
