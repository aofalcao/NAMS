"""
The NAMS python package calculates the similarity between molecules based on the 
structural/topological relationships of each atom towards all the others 
within a molecule.

This program is free software: you can redistribute it and/or modify
it under the terms of the MIT License as published on the official site of Open Source Initiative
and attached above.

Copyright (C) 2013, Andre Falcao and Ana Teixeira, University of Lisbon - LaSIGE

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

Please cite the authors in any work or product based on this material:

AL Teixeira, AO Falcao. 2013. A non-contiguous atom matching structural similarity function. J. Chem. Inf. Model. DOI: 10.1021/ci400324u.

"""


try:
    reload(recoder)
except:
    import recoder
    
try:
    reload(pybel)
except:
    import pybel 

try:
    reload(openbabel)
except:
    import openbabel

def mol2svg(mname, mol ):
    conv = openbabel.OBConversion()
    conv.SetOutFormat("svg")
    conv.AddOption("d", conv.OUTOPTIONS, "")
    conv.AddOption("r", conv.OUTOPTIONS, "1")
    conv.AddOption("c", conv.OUTOPTIONS, "1")
    #this one will add the atom indices
    conv.AddOption("i", conv.OUTOPTIONS, "1")
    #this one will print it in black and white
    conv.AddOption("u", conv.OUTOPTIONS, "1")

    conv.WriteFile(mol.OBMol, "./"+mname+".svg")
    conv.CloseOutFile()


def recode_file(fname_in, fname_out, makeSVG=False):
    rec=recoder.Recoder()
    hs=False
    filin=open(fname_in, "rt")
    lins=filin.readlines()
    filin.close()

    filout=open(fname_out, "wt")
    for lin in lins:
        smi, mol_id = lin.split("\t")
        cid=int(mol_id)
        try:
            mol=pybel.readstring("smi",smi)
            can_smi = mol.write("can").strip()
            if makeSVG==True:
                mol=pybel.readstring("smi",can_smi)
                mol2svg(str(cid), mol)
        except:
            print cid, "oops! cannot convert smiles: ", smiles.strip()
            continue
        hs=True if len(can_smi)==1 else False
        print cid, len(can_smi), mol.molwt, "...",
        #mol, mol_info = rec.get_mol_info("can", can_smi, hs, True)
        mol, mol_info = rec.get_mol_info("can", can_smi, hs, False)  #without chirality
        if mol:
            rec.export_mol_info(mol_info, cid, can_smi, filout, mol.molwt)
        else:
            print "molecule too small", smiles
        print "Done"
    filout.close()
    
if __name__=="__main__":
    import sys
#    try:
        #TODO: all the cmd line processing. Right now it is only reading smiles and integer IDs
    if len(sys.argv)==3:
        recode_file(sys.argv[1], sys.argv[2], makeSVG=False)
    elif len(sys.argv)==4:
        if sys.argv[3]=="-writeimages":
            recode_file(sys.argv[1], sys.argv[2], makeSVG=True)
        else:
            print "No Entiendo"
    else:
        print "No Entiendo"
#    except:
#        print "Usage: makenamsdb mol_file_name nams_file_name [format=[smi|inchi] [colmol=1] [colid=2]"
