
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

class Chirality:
    """
    This class determinates the chirality of a molecule according to the R-S convention based
    on linear notations to represent the molecular structure: SMILES or InCHI.
    It requires Python 2.6 or above and Openbabel 2.3.1.
    """    
    def __init__(self, line_notation, input_type, isotopes=False):
        self.obmol=None
        self.n_atoms=0
        self.all_bonds={}
        self.chiral_atoms=[]
        self.chiralities=[]
        
        
        obConversion = openbabel.OBConversion()
        mol_can = openbabel.OBMol()
        obConversion.SetInAndOutFormats(input_type, "can")
        obConversion.ReadString(mol_can, line_notation)
        can_smi = obConversion.WriteString(mol_can)
        self.can_smi = can_smi
        #set obmol with Hs based on the canonical smiles
        obConversion.SetInFormat("can")
        mol = openbabel.OBMol()
        obConversion.ReadString(mol, can_smi)
        mol.AddHydrogens()
        self.obmol=mol
        #How many chiral atoms (@) are referenced in the canonical smile
        #and append the result to chiral_atoms
        chir=False
        for c in can_smi:
            if c=="@":
                if chir==True:
                    self.chiral_atoms.append(1)
                    chir=False
                else:
                    chir=True
            else:
                if chir==True:
                    self.chiral_atoms.append(-1)
                    chir=False
        
        self.n_atoms=self.obmol.NumAtoms()
        self.n_bonds = self.obmol.NumBonds()
        self.all_bonds={}
        self.auto_chirs=0
        for atom_id in range(self.n_atoms):
            at=self.obmol.GetAtom(atom_id+1)
            #check how many atoms are INDEED chiral
            if at.IsChiral(): self.auto_chirs+=1
            
            bonds = self.get_bonds(self.obmol, atom_id)
            for b in bonds:
                self.all_bonds.setdefault(b[0], {})[b[1]]=1
                self.all_bonds.setdefault(b[1], {})[b[0]]=1
        
        #test the case: the molecule is chiral but it is not completely
        #stated, so it is ignored
        #1 corresponds to R and -1 corresponds to S
        if len(self.chiral_atoms)==self.auto_chirs: self.calc_chiralities(isotopes)
        else:
            ##verify stereo bonds
            self.chiralities=[0]*self.n_atoms
        #print "Chirality"
        #print self.chiralities
        
        
        
    def get_chirality(self, atom_id):
        #atom_id is the index of the atom in the molecule (1:nr of atoms)
        #we need the vector index, so it is the atom id of the molecule - 1
        return self.chiralities[atom_id] 
    

    def calc_chiralities(self, isotopes):
        #1. check the chiral atoms - Q
        #2. For each substituent of each chiral atom remove the bonds to Q.
        #3. Process the bonds recursively
        #4. store, in decreasing order, the tuples
        chir_atom=0 
        for atom_id in range(self.n_atoms):
            at=self.obmol.GetAtom(atom_id+1)
            if at.IsChiral():
                #check all neighbours
                chiral_levels=[]
                level0=[]    #this will store the atomic numbers of the direct ligands to the chiral atom
                for n_atom in openbabel.OBAtomAtomIter(at):
                    id2 = n_atom.GetIndex()
                    #add them to the neighbours list
                    elliminated=[]
                    elliminated.append((atom_id, id2))
                    elliminated.append((id2, atom_id))    
                    #get all the bonds for each separation level
                    lvl_bonds=[[]]
                    lvl=0
                    start_set=[id2]
                    self.process_bonds(start_set, elliminated, lvl_bonds)
                    lvl_bonds=lvl_bonds[:-1]
                    at2=self.obmol.GetAtom(id2+1)
                    #level0 is the atomic number of the ligand
                    #the other levels are the level_bonds
                    levels=[]
                    if isotopes:
                        levels.append(at2.GetExactMass())
                    else:
                        levels.append(at2.GetAtomicNum())
                    for lvl in lvl_bonds:
                        levels.append(lvl)
                    chiral_levels.append(levels)
                chiral=self.is_chiral(chiral_levels,  chir_atom)
                self.chiralities.append(chiral)
                chir_atom += 1
            else:
                self.chiralities.append(0)


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
        #for a given atom in a molecule finds the bonds
        natoms=mol.NumAtoms()
        bonds=[]
        obatom = mol.GetAtom(atom_id+1)
        id1 = obatom.GetIndex()
        for n_atom in openbabel.OBAtomAtomIter(obatom):
            id2 = n_atom.GetIndex()
            bonds.append((id1, id2))
        return bonds


    def is_chiral(self, chiral_levels, chir_atom):
        #it will return for the atom in cause whether it is R or S
        #irotation will indicate the sense of rotation of the initial order
        irotation= self.chiral_atoms[chir_atom]
        magic_list=[]
        i=0
        for lig in chiral_levels:
            lig_id=[[lig[0]]]
            for lvl in lig[1:]:
                nlvl=[]
                for bnd in lvl:
                    #will get the second atom in the bond and get the
                    #respective atomic number
                    at=self.obmol.GetAtom(bnd[1]+1)
                    nlvl.append(at.GetAtomicNum())
                nlvl.sort()
                nlvl.reverse()
                lig_id.append(nlvl)
            magic_list.append((lig_id, i))
            i+=1
        #order list according to the CIP rules
        magic_list.sort()
        permutations = [] 
        duplicate_ligs={}
        aux=[x[0] for x in magic_list]
        for lig in magic_list:
            if aux.count(lig[0]) > 1:
                duplicate_ligs[lig[1]] = chiral_levels[lig[1]] #lig[0]           
            permutations.append(lig[1])
        #order the permutation matrix
        permutations.reverse()
##      CIP Rule nr 5 - not used in the current version (0.9)        
        if duplicate_ligs:
            return 0
        return self.parity(permutations, irotation)

    
    def priority_e_z(self, duplicate_ligs, permutations):

        #test for double bonds priority - rule 5
        #it is not being used in the current version (0.9)
        
        import doubleb_e_z
        ez=doubleb_e_z.Stereodoubleb(self.can_smi)
        e_z_lig={}
        print duplicate_ligs
        for k in duplicate_ligs:             
            for x in duplicate_ligs[k][1:len(duplicate_ligs[k])]:
               for (at1, at2) in x:
                   at_1=self.obmol.GetAtom(at1+1)
                   at_2=self.obmol.GetAtom(at2+1)
                   bond = at_2.GetBond(at_1)
                   if bond.IsDouble():
                       b_idx = bond.GetIdx()
                       if ez.get_e_z(at_1, at_2) != 0:
                           if k in e_z_lig.keys():  
                               e_z_lig[k].append(ez.get_e_z(at_1, at_2))  
                           else:  
                               e_z_lig[k] = [ez.get_e_z(at_1, at_2)]                  
        keys = e_z_lig.keys()
        import operator
        sorted_x = sorted(e_z_lig.iteritems(), key=operator.itemgetter(1))
        x = 0
        permutations_old = permutations[:]
        for ezlig in sorted_x:
            idx = permutations_old.index(keys[x])
            permutations[idx] = ezlig[0]
            x = x + 1
            
        return permutations


    def priority_chiral(self, duplicate_ligs, permutations):
        #test for chiral centers in the duplicate ligands - rule 5
        #it is not being used in the current version (1.0)
        chir_lig={}
        import chiraux
        #ch = chiraux.Chirality(self.obmol, self.can_smi)
        for k in duplicate_ligs:
            #find unique atoms
            #find chiral
            #get_chirality
            unique_atoms=[]
            for x in duplicate_ligs[k][1:len(duplicate_ligs[k])]:
               for (at1, at2) in x:
                   if at1 not in unique_atoms:
                       unique_atoms.append(at1)
                       at_1=self.obmol.GetAtom(at1+1)
                       if at_1.IsChiral():
                           if k in chir_lig.keys():
                               chir_lig[k].append(ch.get_chirality(at1))
                               #chir_lig[k].append(self.get_chirality(at1))
                           else:
                               chir_lig[k] = [ch.get_chirality(at1)]
                               #chir_lig[k] = [self.get_chirality(at1)] 
                   if at2 not in unique_atoms:
                       unique_atoms.append(at2)
                       at_2=self.obmol.GetAtom(at2+1)
                       if at_2.IsChiral():
                           if k in chir_lig.keys():  
                               chir_lig[k].append(ch.get_chirality(at2))
                               #chir_lig[k].append(self.get_chirality(at2))
                           else:  
                               chir_lig[k] = [ch.get_chirality(at2)]
                               #chir_lig[k] = [self.get_chirality(at2)] 
        keys = chir_lig.keys()
        import operator
        sorted_x = sorted(chir_lig.iteritems(), key=operator.itemgetter(1))
        x = 0
        permutations_old = permutations[:]
        for chirlig in sorted_x:
            idx = permutations_old.index(keys[x])
            permutations[idx] = chirlig[0]
            x = x + 1
        return permutations
        
    def parity(self, permutation, irotation):
        #in order to determine whether a given permutation is even or odd
        #one writes the permutation as a product of disjoint cycles.
        #The permutation is odd if and only if this factorization contains an *odd number* of *even-length* cycles.
        #cycles with 1 element do not count.
        # a cycle ends when the 1st element of the cycle is found
        #example from [0 1 2 3] permutation = [1 3 2 0]: 0->1;1->3;2->2;3->0 => (013)(2) -> odd number (1) of odd-length(3) cycles ->parity = EVEN!
        #example from [0 1 2 3] permutation = [1 2 3 0]: 0->1;1->2;2->3;3->0 => (0123) -> odd number (1) of even-length (4) cycles ->parity = ODD!
        i = 1
        e = permutation[0]
        #e = i!
        #Find the position of "0" in the permutation in order to verify the length
        #of the first cycle (normally it ends with zero, except when zero is in the 1st position)
        #the position of 0 indicates us the size of the first cycle 
        while e != 0:
           #finds the next element of the cycle
           e = permutation[e]
           i += 1

        #one cycle with length 2, we have 2 possibilities for the remaining: 2 cycles with length 1 or 1 cycle with length 2.
        #possibilities:(X0)(X)(X)->odd number(1)of even length cycles(2)=> parity = ODD!
        #(X0)(XX) -> even number (2) of even-length cycles (2)=>parity = EVEN!
        if i == 2:
           # permutation[0] is  1, 2, or 3.
           # we need to know if the second cycle has 2 or 1 elements
           #aux will save the first element of the first cycle (we are not interested)
           
           aux = permutation[0]
           zero = permutation[aux]
           #we should now test if the element is in the correct position, ie if the element in the permutation is in the same
           #position of the original order [0 1 2 3 4] -> this means we have 2 cycles of length one, else it is 1 cycle of length 2.
           #lets search for the other elements
           for i in range(0,3):
               if permutation[i] != aux and permutation[i] != zero:
                   #then we have 2 cycles of lenght 1 => parity=ODD!
                   if permutation[i] == i and irotation == 1:
                       return -1 # "odd"
                   elif permutation[i] == i and irotation == -1:
                       return 1      # "even"
                   elif permutation[i] != i and irotation == 1:
                       return 1 #even
                   elif permutation[i] != i and irotation == -1:
                       return -1 #odd

        #One cycle of length 3 means the other is length 1
        #zero is in the 3rd position -> this always generates 1 cycles (odd) with length 3 (XX0)(X) (odd) => parity=EVEN! 
        if i == 3:
            #if irotation is @ then it is even
            if irotation == 1:
               return 1 #"even"
            #else it is odd
            elif irotation == -1:
                return -1 #"odd"
            
        #Everything is involved in a cycle of length 4
        #zero is in the 4th position -> this always generates 1 cycle (odd) with length 4 (XXX0) (even) => parity=ODD! 
        if i == 4:
           #if irotation is @ then it is odd
           if irotation == 1:
               return -1 #"odd"
            #else it is even
           elif irotation == -1:
               return 1 #"even"

        #zero is in the first position, so we should start searching for the element 1 and from the position 1. 
        # permutation[0] == 0, so check the cycle starting at permutation[1]
        i = 1
        e = permutation[1]
        while e != 1:
           e = permutation[e]
           i += 1

        #the cycle has length 2 -> this means we have (X)(X1)(X) -> this generates 1 cycle (odd) with length 2(even)=>parity=ODD!
        if i == 2:
           #if irotation is @ then it is odd
           if irotation == 1:
               return -1 #"odd"
           #else it is even
           elif irotation == -1:
               return 1  #"even"
            
        #the cycle has length 3 -> this means we have (X)(XX1)-> this generates 1 cycle (odd) with length 3(odd)=>parity=EVEN!   
        if i == 3:
           #if irotation is @ then it is even
           if irotation == 1:
               return 1  #"even"
           #else it is odd
           elif irotation == -1:
               return -1  #"odd"
        #the last case is: zero in the 1st position and 1 in the 2nd position
        #we have 2 possibilities: 2 cycles of length 1 or 1 cycle of length 2.
        # then we will check if the element in the 3rd position of the array (2) corresponds to 2, if it does then we have 2 cycles of
        # length 1
        # permutation[0] == 0, permutation[1] == 1, what about permutation[2]?
        if permutation[2] == 2 and irotation == 1:
            # this means we have (1)(2)(3)(4) - This is the identity =>parity=EVEN!
           return 1 #"even"
        elif permutation[2] == 2 and irotation == -1:
            return -1 #odd
        #this means we have (1)(2)(34) -> 1 cycle (odd) with length 2(even)=>parity=ODD!
        elif permutation[2] != 2 and irotation == 1:
            return -1   #"odd"            
        elif permutation[2] != 2 and irotation == -1:
            return 1   #"even"

