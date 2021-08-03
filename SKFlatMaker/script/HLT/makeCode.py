import argparse

parser=argparse.ArgumentParser()
parser.add_argument("file")
args=parser.parse_args()

filters=[]
with open(args.file) as f:
    for line in f.readlines():
        line=line.split("#")[0]
        line=line.strip()
        if line=="": continue
        filters+=[line]
filters=sorted(set(filters))

print "for(pat::TriggerObjectStandAlone obj : *triggerObject){"
print "  if(deltaR(obj,imuon)<0.3){"
print "    obj.unpackFilterLabels(iEvent, *trigResult);"
for i in range(len(filters)):
    print "    if(obj.filter(\"{}\")) filterbits|=1<<{};".format(filters[i],i)
print "  }"
print "}"
