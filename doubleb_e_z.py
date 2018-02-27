

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

import openbabel, pybel
from math import atan, sin, cos

class Stereodoubleb:
    
    def __init__(self, line_notation, input_type, convtype="OB", isotopes=False):
        self.obmol=None
        self.n_bonds=0
        self.n_atoms=0
        self.all_bonds={}
        self.double_b=[]
        self.e_z=[]
        self.isotopes = isotopes
        self.number_explicit_db = 0
        obConversion = openbabel.OBConversion()
        mol_can = openbabel.OBMol()
        obConversion.SetInAndOutFormats(input_type, "can")
        obConversion.ReadString(mol_can, line_notation)
        can_smi = obConversion.WriteString(mol_can)
        self.can_smi = can_smi
        
        #sets the molecule - mol MUST BE A OBMOL
        #using the canonical SMILES
        obConversion = openbabel.OBConversion()
        obConversion.SetInFormat("can")
        mol = openbabel.OBMol()
        obConversion.ReadString(mol, can_smi)
        mol.AddHydrogens()
        self.obmol=mol
        self.n_bonds=self.obmol.NumBonds()
        self.all_bonds={}
        double_b=0

        
        for bond_id in range(self.n_bonds):
            bond=self.obmol.GetBond(bond_id)
            #check how many double bonds exist in the molecule
            if bond.IsDouble() and not bond.IsAromatic():
                double_b+=1
                
        if double_b>0:
            if convtype=="MARVIN":
                #use Marvin - molconvert (ChemAxon) to produce the MOLFILE
                #WARNING: indexes may be different
                molfile = convert_molconv("can", can_smi)
                #sets the molecule with Hs and 2D coordinates
                obConversion = openbabel.OBConversion()
                mol = openbabel.OBMol()
                obConversion.SetInFormat("mol")
                obConversion.ReadString(mol, molfile)
                self.obmol=mol
            else:
                #use openbabel to produce the 2D MOL
                self.obmol=convert(mol)

            self.number_explicit_db = explicit_cistrans(self.can_smi)
            self.n_bonds=self.obmol.NumBonds()
            self.n_atoms=self.obmol.NumAtoms()
            for atom_id in range(self.n_atoms):
                at=self.obmol.GetAtom(atom_id+1) 
                bonds = self.get_bonds(self.obmol, atom_id)
                for b in bonds:
                    self.all_bonds.setdefault(b[0], {})[b[1]]=1
                    self.all_bonds.setdefault(b[1], {})[b[0]]=1
            self.calc_e_z()
        else:
            self.e_z=[0]*self.n_bonds
        implicit_db = 0
        for x in self.e_z:
            if x != 0:
                implicit_db = implicit_db + 1
        if self.number_explicit_db != implicit_db:
            self.e_z=[0]*self.n_bonds

    def get_e_z_at(self, at1, at2):
        #receive the atoms
        #bond_id is the bond id of the molecule
        bond_id = at2.GetBond(at1).GetIdx()
        return self.e_z[bond_id]
    
    def get_e_z_bond(self, bond_idx):
        #receive the bond index
        return self.e_z[bond_idx]  


    def calc_e_z(self):
        #1. check the double bonds and the 2 atoms forming it
        #2. save the coordinates of the 2 atoms forming the double bonds
        #3. For each substituent of each carbon atom of the double bond remove the bonds to it.
        #3. Process the bonds recursively
        #4. store, in decreasing order, the tuples       
        double_b=0 
        for bond_id in range(self.n_bonds):
            bond=self.obmol.GetBond(bond_id)
            at1 = bond.GetBeginAtom()
            at2 = bond.GetEndAtom()
            if bond.IsDouble() and not bond.IsAromatic() and at1.IsCarbon() and at2.IsCarbon() and at1.CountBondsOfOrder(2) == 1 and at2.CountBondsOfOrder(2)== 1:
                #check all neighbours of at1 and at2 (at1=at2)
                coordinatesc =[]
                coordinatesc.append([at1.GetX(), at1.GetY()])
                coordinatesc.append([at2.GetX(), at2.GetY()])
                self.coordinatesc = coordinatesc
                level0=[]    #this will store the atomic numbers of the direct ligands to
                #each atom of the double bond atom
                atoms_indoubleb =[at1, at2]
                coordinates = []
                double_atom = []
                #process the substituents of each atom of the double bond
                
                x = 0
                for x in range (len(atoms_indoubleb)):
                    at =  atoms_indoubleb[x]
                    if x == 0:
                        atom_2 = atoms_indoubleb[x+1]
                    else:
                        atom_2 = atoms_indoubleb[x-1]
                    coordinates.append([])
                    double_atom.append([])
                    for n_atom in openbabel.OBAtomAtomIter(at):
                        id2 = n_atom.GetIndex()
                        if id2 != atom_2.GetIndex():
                            #add them to the neighbours list
                            elliminated=[]
                            elliminated.append((at.GetIndex(), id2))
                            elliminated.append((id2, at.GetIndex()))    
                            #get all the bonds for each separation level
                            lvl_bonds=[[]]
                            lvl=0
                            start_set=[id2]
                            self.process_bonds(start_set, elliminated, lvl_bonds)
                            lvl_bonds=lvl_bonds[:-1]
                            atom2=self.obmol.GetAtom(id2+1)
                            #level0 is the atomic number of the ligand
                            #the other levels are the level_bonds
                            levels=[]
                            if self.isotopes:
                                levels.append(atom2.GetExactMass())
                            else:
                                levels.append(atom2.GetAtomicNum())
                            for lvl in lvl_bonds:
                                levels.append(lvl)
                            double_atom[x].append(levels)
                            coordinates[x].append([n_atom.GetX(), n_atom.GetY()])
                    x = x + 1
                self.coordinates = coordinates
                self.double_atom = double_atom
                verify_e_z=self.is_e_z()
                self.e_z.append(verify_e_z)
                double_b += 1
            else:
                self.e_z.append(0)


    def process_bonds(self, start_set, elliminated, lvl_bonds):
        #the idea is to select all the immediate bonds in the start set
        #append them to the current level_bonds and elliminate them from the
        #connection matrix
        to_follow={}
        for at1 in start_set:
            for at2 in self.all_bonds[at1]:
                if (at1,at2) not in elliminated:
                    my_bond=(at1, at2)
                    b=self.obmol.GetBond(at1+1, at2+1)
                    bo=b.GetBO()
                    for i in range(0, bo):
                        lvl_bonds[-1].append(my_bond)
                        
                    if at2<self.n_atoms: to_follow.setdefault(at2,1)
                    elliminated.append((at1,at2))
                    elliminated.append((at2,at1))
        #pick up the new nodes (from the bonds) and create a new start set
        #process recursively until no more nodes are extant to process
        #the recursion proceeds at most as many times as the total length
        #of the graph 
        new_starts=to_follow.keys()
        if len(new_starts)>0:
            lvl_bonds.append([])
            self.process_bonds(new_starts, elliminated, lvl_bonds)


    def get_bonds(self, mol, atom_id):
        #for a given atom in a molecule find the bonds
        natoms=mol.NumAtoms()
        bonds=[]
        obatom = mol.GetAtom(atom_id+1)
        id1 = obatom.GetIndex()
        for n_atom in openbabel.OBAtomAtomIter(obatom):
            id2 = n_atom.GetIndex()
            bonds.append((id1, id2))
        return bonds


    def is_e_z(self):
        #here it will return for the bond in cause whether it is E or Z
        #the rotated coordinates will indicate if the higher priority group
        #on each side of the double bond is up or down 
        #Up = 1
        #Down = -1
        #0 - groups are equal
        
        vertical = False #indicates if the molecule is in vertical position or not
        if self.coordinatesc[0][0]-self.coordinatesc[1][0] != 0:
            #calculate the slope of the double bond, based on the x,y coordinates of the C atoms
            slope = (self.coordinatesc[1][1] - self.coordinatesc[0][1])/(self.coordinatesc[1][0] - self.coordinatesc[0][0])
            #calculate the rotation angle based on the following expression:
            #teta=tan-1[(yo-y1)/(x0-x1)]
            teta=atan(slope)
            #if slope is positive, rotate the molecule clockwise
            if slope > 0:
                teta=-teta              
        else:
            vertical = True
        u_d=[]
        

        for d_a in range(0, len(self.double_atom)):
            da_ligs = self.double_atom[d_a]
            magic_list=[]
            i = 0
            for lig in da_ligs:
                lig_id=[[lig[0]]]       
                for lvl in lig[1:]:
                    nlvl=[]
                    for bnd in lvl:
                        #will get the second atom in the bond and get the respective atomic number
                        at=self.obmol.GetAtom(bnd[1]+1)
                        nlvl.append(at.GetAtomicNum())
                    nlvl.sort()
                    nlvl.reverse() 
                    lig_id.append(nlvl)
                magic_list.append((lig_id, i))
                i+=1
            
            magic_list.sort()
            magic_list.reverse()
            if magic_list[0][0] == magic_list[1][0]:
                return 0
            
            #x=x*sin(teta)+ycos(teta)
            if vertical:
                #compare X
                if self.coordinates[d_a][magic_list[0][1]][0] > self.coordinates[d_a][magic_list[1][1]][0]:
                    u_d.append(1)
                else:
                    u_d.append(-1) 
            else:
               
                y1=self.coordinates[d_a][magic_list[0][1]][0]*sin(teta)+self.coordinates[d_a][magic_list[0][1]][1]*cos(teta)
                y2=self.coordinates[d_a][magic_list[1][1]][0]*sin(teta)+self.coordinates[d_a][magic_list[1][1]][1]*cos(teta)

                if y1 > y2:
                    u_d.append(1)
                else:
                    u_d.append(-1)


        #Possibilities: cis: (1)(1); (-1)(-1); trans: (1)(-1); (-1)(1)
        if u_d[0] == u_d[1]:
            return 1 #cis - Z
        else:
            return -1 #trans - E
        
def explicit_cistrans(can_smi):
    obConversion = openbabel.OBConversion()
    mol = openbabel.OBMol()
    obConversion.SetInAndOutFormats("can", "inchi")
    obConversion.ReadString(mol, can_smi)
    inchi = obConversion.WriteString(mol)
    recon_start = inchi.find("/b")
    doub_info=[]
    if recon_start != -1:
        doub_info = inchi[recon_start+2:].rstrip().replace(";",",").split(",")

    return len(doub_info)
  
                    
def convert(molecule):
    gen2d = openbabel.OBOp.FindType("Gen2D")
    mol=molecule
    gen2d.Do(molecule)
    return molecule
    

def convert_molconv(tipo, molecule):
    import os
    path = "C:\\Program Files (x86)\\ChemAxon\\MarvinBeans\\bin"
    torun = "molconvert mol:H -s \"" + molecule.rstrip() + "\""
    os.chdir(path)
    w,r=os.popen2(torun)
    molf=r.readlines()
    molfile=""
    for i in molf:
        molfile += i
    return molfile
