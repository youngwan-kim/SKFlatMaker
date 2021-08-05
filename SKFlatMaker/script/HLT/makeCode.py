import argparse

parser=argparse.ArgumentParser()
args=parser.parse_args()

def readFile(filename):
    filters=[]
    for line in open(filename).readlines():
        line=line.split("#")[0]
        line=line.strip()
        if line=="": continue
        filters+=[line]
    filters=sorted(set(filters))
    return filters
    
muon_filters=readFile("muon_filters.txt")
muon_paths=readFile("muon_paths.txt")
print "------------------ Muon -----------------"
print "for(pat::TriggerObjectStandAlone obj : *triggerObject){"
print "  if(deltaR(obj,imuon)<0.3){"
print "    obj.unpackPathNames(trigNames);"
print "    obj.unpackFilterLabels(iEvent, *trigResult);"
for i in range(len(muon_paths)):
    print "    if(obj.path(\"{}*\",true,true)) pathbits|=ULong64_t(1)<<{};".format(muon_paths[i],i)
for i in range(len(muon_filters)):
    print "    if(obj.filter(\"{}\")) filterbits|=ULong64_t(1)<<{};".format(muon_filters[i],i)
print "  }"
print "}"

electron_filters=readFile("electron_filters.txt")
electron_paths=readFile("electron_paths.txt")
print "------------------ Electron -----------------"
print "for(pat::TriggerObjectStandAlone obj : *triggerObject){"
print "  if(deltaR(obj,*el)<0.3){"
print "    obj.unpackPathNames(trigNames);"
print "    obj.unpackFilterLabels(iEvent, *trigResult);"
for i in range(len(electron_paths)):
    print "    if(obj.path(\"{}*\",true,true)) pathbits|=ULong64_t(1)<<{};".format(electron_paths[i],i)
for i in range(len(electron_filters)):
    print "    if(obj.filter(\"{}\")) filterbits|=ULong64_t(1)<<{};".format(electron_filters[i],i)
print "  }"
print "}"
