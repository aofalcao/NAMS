
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
import pybel, openbabel
import sys
from chirality import Chirality
from doubleb_e_z import Stereodoubleb
from math import sqrt
import os


class Recoder:

    def get_num_rings(self, mol, atom):
        #counts the number of rings an atom belongs to
        cont_rings = 0
        for ring in mol.GetSSSR():
            if ring.IsMember(atom):
                cont_rings += 1
        
        return cont_rings



    def process_bonds(self, start_set, all_bonds, n_big_atoms, elliminated, lvl_bonds):
        #the idea is to select all the immediate bonds in the start set
        # append them to the current level_bonds and elliminate them from the connection matrix
        to_follow={}
        for at1 in start_set:
            for at2 in all_bonds[at1]:
                if (at1, at2) not in elliminated:
                    my_bond=(at1, at2)
                                  
                    lvl_bonds[-1].append(my_bond)
                    
                    if at2<=n_big_atoms: to_follow.setdefault(at2,1)
                    elliminated.append((at1,at2))
                    elliminated.append((at2,at1))
        # now pick up the new nodes (from the bonds) and create a new start set
        #process recursively until no more nodes are extant to process
        #the recursion proceeds at most as many times as the total length of the graph 
        new_starts=to_follow.keys()

        if len(new_starts)>0:
            lvl_bonds.append([])
            self.process_bonds(new_starts, all_bonds, n_big_atoms, elliminated, lvl_bonds)
                


    def calc_btypes(self, mol_info):
        # get the existing bondtypes from a mol_info
        # output a dictionary with the bond types
        #     and the mol info redesigned for the bond types which can then be forgotten
        btypes={}
        nbonds={}
        idx=0
        for at in mol_info:
            lvls = mol_info[at]
            nbonds[at]=[]
            L=0
            for lvl in lvls:
                for bond in lvl:
                    if bond not in btypes:
                        btypes[bond]=idx
                        nbonds[at].append((idx, L))
                        idx+=1
                    else:
                        nbonds[at].append((btypes[bond], L))
                L+=1
        return btypes, nbonds



    def get_bonds2(self, mol, atom_id):
        #for a given atom in a molecule find the bonds
        natoms=mol.NumAtoms()
        bonds=[]
        obatom = mol.GetAtom(atom_id)
        id1 = obatom.GetIndex()
        t1 = obatom.GetType()
        for n_atom in openbabel.OBAtomAtomIter(obatom):
            id2 = n_atom.GetIndex()
            t2=n_atom.GetType()
            bond = obatom.GetBond(n_atom)
            bonds.append((id1+1, id2+1))
        return bonds


    def get_mol_info(self, typ, mol_str, hydros=False, DoIsomerism=False):
        
        mol = pybel.readstring(typ, mol_str)
        obmol=mol.OBMol
        n_big_atoms = len(mol.atoms)
        if n_big_atoms<2: return False, False
		
        #use canonical smiles ONLY!!!
        can_smi = mol.write("can")
        #remove other atoms not bonded to the molecule
        if can_smi.find("."): can_smi=can_smi.split(".")[0]

        mol = pybel.readstring("smi", can_smi)
        obmol=mol.OBMol
        n_big_atoms = len(mol.atoms)
        #this is the chirality/double bond stereo stuff (both work with explicit Hs)
        if DoIsomerism==True:
            chir=Chirality(can_smi, "can")
            dbstereo=Stereodoubleb(can_smi, "can")
               
        if hydros==True: obmol.AddHydrogens()
        n_tot_atoms=obmol.NumAtoms()
        if n_big_atoms<2:
            return False, False
	all_bonds={}
        mol_info={}
        mol_atoms=[]
        

        #start up the connection matrix as a dictionary of dictionaries (sparse matrix)
        
        for atom_id in range(n_tot_atoms):
            ida=atom_id+1
            atm=obmol.GetAtom(ida)
            if atom_id<n_big_atoms and DoIsomerism==True:
                the_chir=chir.get_chirality(atom_id)
            else:
                the_chir=0
            mol_atoms.append((atm, atm.GetAtomicNum(), self.get_num_rings(obmol,atm), the_chir))
            bonds = self.get_bonds2(obmol, atom_id+1)
            for b in bonds:
                all_bonds.setdefault(b[0], {})[b[1]]=1   #b[2]
                all_bonds.setdefault(b[1], {})[b[0]]=1   #b[2]
        #print n_big_atoms
        #for ab in all_bonds: print ab, all_bonds[ab]
        #for each BIG atom in the molecule compute the bonds according to their separation levels 
        for atom_id in range(1,n_big_atoms+1):
            elliminated=[]
            lvl_bonds=[[]]
            lvl=0
            start_set=[atom_id]
            #get all the bonds for each separation level
            self.process_bonds(start_set, all_bonds, n_tot_atoms, elliminated, lvl_bonds)
            #for some weird reason the last level is always empty so take it out
            #for lvl in lvl_bonds: print lvl
            lvl_bonds=lvl_bonds[:-1]
            the_bonds=[]
            #print "-->", atom_id
            for level in lvl_bonds:
                #print "\t-->", level
                the_bonds.append([])
                for at1, at2 in level:
                    #print "\t\t-->", at1, at2
                    ai1, ai2=at1-1, at2-1
                    a1 = mol_atoms[ai1][0]
                    a2 = mol_atoms[ai2][0]
                    
                    bnd12 = a1.GetBond(a2)
                    #get double bond stereoisomerism
                    dbstereo12=0
                    if DoIsomerism==True:
                        dbstereo12 = dbstereo.get_e_z_at(a1, a2)
                    if bnd12.IsAromatic():
                        my_bond=(mol_atoms[ai1][1], mol_atoms[ai1][2], mol_atoms[ai1][3],
                                 mol_atoms[ai2][1], mol_atoms[ai2][2], mol_atoms[ai2][3],
                                 bnd12.IsInRing(),bnd12.IsAromatic(), 1.5, dbstereo12) 
                    else:
                        my_bond=(mol_atoms[ai1][1], mol_atoms[ai1][2], mol_atoms[ai1][3],
                                 mol_atoms[ai2][1], mol_atoms[ai2][2], mol_atoms[ai2][3],
                                 bnd12.IsInRing(),False, bnd12.GetBondOrder(), dbstereo12)  

                    
                    the_bonds[-1].append(my_bond)
                    
            mol_info[atom_id]=the_bonds
        return mol, mol_info

    def export_mol_info(self, mol_info, cid, cname, fil, molwt=-1.0):
        ababs={}
        #this function will export in a text file the mol_info requiired for nams
        #print mol_info
        natoms = len(mol_info)
        #here we will get all the different aba_bonds
        new_id=0
        mat_aba_typs=[]
        mat_levels=[]
        for a in mol_info.keys():
            mat_aba_typs.append([])
            mat_levels.append([])
            nlevels=len(mol_info[a])
            #print a, nlevels
            nbonds=0
            for L in range(nlevels):
                level=mol_info[a][L]
                nbonds+=len(level)
                for abab in level:
                    #print "\t", L, abab
                    if abab not in ababs:
                        ababs[abab]=new_id
                        new_id+=1
                    aba_typ=ababs[abab]
                    mat_aba_typs[-1].append(aba_typ)
                    mat_levels[-1].append(L)
        fil.write("%d %d %s\n" % (cid, int(molwt*10), cname))
        fil.write("%d %d %d\n" %(natoms, nbonds, len(ababs)))
        #print them
        #sort them according to value using a zip - this is just for output!
        Z=zip(ababs.values(), ababs.keys())
        Z.sort()
        for v, aba in Z:
            #print aba, v
            #float numb1, nrings1, chir1, numb2, nrings2, chir12, inring, arom, order, dbcistrans;
            fil.write("%d %d %d %d %d %d %d %d %d %d\n" %(aba[0], aba[1], aba[2], aba[3], aba[4],
                                                    aba[5], aba[6], aba[7], aba[8]*10, aba[9]))
        #now process the matrix atoms vs ababonds
        for atom in mat_aba_typs:
            s=reduce(lambda p, q: str(p)+" " +str(q), atom, "")
            fil.write(s.strip()+"\n")
        #finally process the matrix atoms vs levels
        #print
        for atom in mat_levels:
            s=reduce(lambda p, q: str(p)+" " +str(q), atom, "")
            fil.write(s.strip()+"\n")
            
    def export_mol_info_bin(self, mol_info, cid, cname, fil, molwt=-1.0):
        import struct
        ababs={}
        #this function will export in a text file the mol_info requiired for nams
        #print mol_info
        natoms = len(mol_info)
        #here we will get all the different aba_bonds
        new_id=0
        mat_aba_typs=[]
        mat_levels=[]
        for a in mol_info.keys():
            mat_aba_typs.append([])
            mat_levels.append([])
            nlevels=len(mol_info[a])
            #print a, nlevels
            nbonds=0
            for L in range(nlevels):
                level=mol_info[a][L]
                nbonds+=len(level)
                for abab in level:
                    #print "\t", L, abab
                    if abab not in ababs:
                        ababs[abab]=new_id
                        new_id+=1
                    aba_typ=ababs[abab]
                    mat_aba_typs[-1].append(aba_typ)
                    mat_levels[-1].append(L)
        #fil.write("%d %s\n" % (cid, cname))
        #fil.write("%d %d %d\n" %(natoms, nbonds, len(ababs)))
        #the molecular weight is FINALLY included
        fil.write(struct.pack("II32s", cid, int(molwt*10), cname))
        fil.write(struct.pack("hhh", natoms, nbonds, len(ababs)))
        #print them
        #sort them according to value using a zip - this is just for output!
        Z=zip(ababs.values(), ababs.keys())
        Z.sort()
        for v, aba in Z:
            #print aba, v
            #float numb1, nrings1, chir1, numb2, nrings2, chir12, inring, arom, order, dbcistrans;
            #fil.write("%d %d %d %d %d %d %d %d %d %d\n" %(aba[0], aba[1], aba[2], aba[3], aba[4],
            #                                              aba[5], aba[6], aba[7], aba[8]*10, aba[9]))
            fil.write(struct.pack("BBBBBBBBBB", aba[0], aba[1], aba[2], aba[3], 
                                  aba[4], aba[5], aba[6], aba[7], aba[8]*10, 
                                  aba[9]))
        #now process the matrix atoms vs ababonds
        for atom in mat_aba_typs:
            for i in atom:
                fil.write(struct.pack("B", i))
                
            #s=reduce(lambda p, q: str(p)+" " +str(q), atom, "")
            #fil.write(s.strip()+"\n")
        #finally process the matrix atoms vs levels
        #print
        for atom in mat_levels:
            for i in atom:
                fil.write(struct.pack("B", i))
            #s=reduce(lambda p, q: str(p)+" " +str(q), atom, "")
            #fil.write(s.strip()+"\n")
        

            
        

