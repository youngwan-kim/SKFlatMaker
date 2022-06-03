import sys
import ROOT

f1=ROOT.TFile(sys.argv[1])
f2=ROOT.TFile(sys.argv[2])
tree1=f1.Get("recoTree/SKFlat")
tree2=f2.Get("recoTree/SKFlat")
simpletypes1=['Int_t', 'ULong64_t', 'Bool_t']
simpletypes2=['Float_t','Double_t']
vectortypes1=['vector<int>', 'vector<string>', 'vector<unsigned int>', 'vector<bool>', 'vector<ULong64_t>']
vectortypes2=['vector<float>','vector<double>']
debug=False
for ie in range(tree1.GetEntries()):
    print(ie)
    tree1.GetEntry(ie)
    tree2.GetEntry(ie)
    for leaf in tree1.GetListOfLeaves():
        name=leaf.GetName()
        tree1attr=getattr(tree1,name)
        tree2attr=getattr(tree2,name)
        if leaf.GetTypeName() in simpletypes1:
            if tree1attr!=tree2attr:
                print ie,name,tree1attr,tree2attr
            if debug: print ie,name,tree1attr,tree2attr
        elif leaf.GetTypeName() in simpletypes2:
            diff=abs(tree1attr-tree2attr)
            if abs(tree2attr): diff/=abs(tree2attr)
            if diff>1e6:
                print ie,name,tree1attr,tree2attr
        elif leaf.GetTypeName() in vectortypes1:
            if tree1attr.size()!=tree2attr.size():
                print ie,name,"size",tree1attr.size(),tree2attr.size()
            for iv in range(tree1attr.size()):
                if type(tree1attr[iv]) is ROOT._Bit_reference: 
                    tree1val=bool(tree1attr[iv])
                else:
                    tree1val=tree1attr[iv]
                if type(tree2attr[iv]) is ROOT._Bit_reference: 
                    tree2val=int(bool(tree2attr[iv]))
                else:
                    tree2val=tree2attr[iv]
                if tree1val!=tree2val:
                    print ie,name,iv,tree1val,tree2val
                if debug: print ie,name,iv,tree1val,tree2val
        elif leaf.GetTypeName() in vectortypes2:
            if tree1attr.size()!=tree2attr.size():
                print ie,name,"size",tree1attr.size(),tree2attr.size()
            for iv in range(tree1attr.size()):
                diff=abs(tree1attr[iv]-tree2attr[iv])
                if abs(tree2attr[iv]): diff/=abs(tree2attr[iv])
                if diff>1e6:
                    print ie,name,iv,tree1attr[iv],tree2attr[iv]
        else:
            print "unknown type",leaf.GetTypeName()

                    
            
