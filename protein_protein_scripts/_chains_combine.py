"""
to combines chains into protein.
Include only term_cutoff modelled residues at the termini.
A modelled loop will not be included if it is longer than loop_cutoff.
If so, "TER" will be inserted to split the chain
"""

from __future__ import print_function

import os
import copy
import glob

MODELLED_PDB_SUBFIX = "_modelled.pdb"
MODELLED_RESIDUES = "REMARK  MODELLED RESIDUES:"
ATOM = "ATOM"
HETATM = "HETATM"
TER = "TER"


class ChainCombine(object):
    def __init__(self, pdb_id, chains, modelling_dir, ions_cofactors_files):
        """
        :param pdb_id: str
        :param chains: list of str; e.g., ["A", "B"]
        :param modelling_dir: str; the directory where modeller was done for missing residues.
        :param ions_cofactors_files: list of str
        """
        self._pdb_id = pdb_id
        self._chains = chains
        self._ions_cofactors_files = ions_cofactors_files

        # dict of dic which is returned  by self._load_chain
        self._original_pdb_data = self._load_chains(pdb_id, chains, modelling_dir)
        
    def trim_residues(self, ter_cutoff, loop_cutoff):
        self._trimed_pdb_data = {}
        for chain in self._chains:
            original_modelled_residues = self._original_pdb_data[chain]["modelled_residues"]
            modelled_residues = self._trim_residues(original_modelled_residues, ter_cutoff, loop_cutoff) 
            self._trimed_pdb_data[chain] = {"modelled_residues":modelled_residues}

            self._trimed_pdb_data[chain]["residues_to_minimize"] = self._residues_to_minimize(modelled_residues)

            self._trimed_pdb_data[chain]["atoms"] = self._trim_atoms(modelled_residues, 
                    self._original_pdb_data[chain]["atoms"], self._original_pdb_data[chain]["residues"])
            self._trimed_pdb_data[chain]["residues"] = self._count_residues(self._trimed_pdb_data[chain]["atoms"])
        return None

    def combine(self):
        atoms = []
        residues_to_minimize = []
        shift = 0
        for chain in self._chains:
            atoms.extend(self._trimed_pdb_data[chain]["atoms"])
            atoms.append("TER")
            residues_to_minimize.extend([ r+shift for r in self._trimed_pdb_data[chain]["residues_to_minimize"] ])
            shift += len(self._trimed_pdb_data[chain]["residues"])

        self._combined_pdb_data = {"atoms":copy.deepcopy(atoms),
                                   "residues_to_minimize":residues_to_minimize, "nresidues":shift}
        self._combined_pdb_data["atoms"] = self._change_resid_sequentially()
        self._combined_pdb_data["natoms"] = self._count_atoms(self._combined_pdb_data["atoms"])

        ions_cofactors = self._load_ions_cofactors(self._ions_cofactors_files)
        if ions_cofactors is not None:
            ic_atoms = []
            ic_natoms = 0
            ic_nresidues = 0
            for chain in ions_cofactors.keys():
                ic_atoms.extend(ions_cofactors[chain]["atoms"])
                ic_natoms += ions_cofactors[chain]["natoms"]
                ic_nresidues += ions_cofactors[chain]["nresidues"]
            
            self._combined_pdb_data["atoms"].extend(ic_atoms)
            self._combined_pdb_data["natoms"] += ic_natoms
            self._combined_pdb_data["nresidues"] += ic_nresidues
        return None

    def write_pdb(self, out=None):
        if out==None:
            out = self._pdb_id + "".join(self._chains) + "_modelled.pdb"

        with open(out, "w") as F:
            header = "REMARK NRESIDUES %d\n"%self._combined_pdb_data["nresidues"]
            header += "REMARK NATOMS %d\n"%self._combined_pdb_data["natoms"]
            header += "REMARK MINIMIZE THESE "
            nres_to_minimize = len(self._combined_pdb_data["residues_to_minimize"])
            for i in range(nres_to_minimize):
                header += "%5d"%self._combined_pdb_data["residues_to_minimize"][i]
                if (i+1)%10 == 0 and i < nres_to_minimize-1:
                    header += "\nREMARK MINIMIZE THESE"
            F.write(header + "\n")
            for line in self._combined_pdb_data["atoms"]:
                F.write(line + "\n")
        return None

    def get_nresidues(self):
        return self._combined_pdb_data["nresidues"]

    def get_natoms(self):
        return self._combined_pdb_data["natoms"]

    def write_residues_to_minimize(self, file):
        nresidues = self._combined_pdb_data["nresidues"]
        out_str = "# nresidues %d\n"%nresidues

        residues = self._combined_pdb_data["residues_to_minimize"]
        residues = ["%d"%r for r in residues]
        out_str += "\n".join(residues)
        open(file, "w").write(out_str)
        return None

    def _load_chains(self, pdb_id, chains, modelling_dir):
        chains_pdb_data = {}
        for chain in chains:
            chains_pdb_data[chain] = self._load_chain(pdb_id, chain, modelling_dir)
        return chains_pdb_data

    def _load_chain(self, pdb_id, chain, modelling_dir):
        """
        :param pdb_id: str
        :param chain: one-letter str
        :param modelling_dir: str, path to modeller results for pdb_id
        :return: dict with keys: "modelled_residues" -> dict {"nter" : [], "loops" : [], "cter" : []}
                                "atoms" -> list of ATOM lines
                                "residues" -> list of residue id
        """
        assert len(chain) == 1, "chain must be a single letter"
        infile = os.path.join(modelling_dir, pdb_id + chain + MODELLED_PDB_SUBFIX)
        pdb_data = {}
        with open(infile, "r") as F:
            for line in F:
                if MODELLED_RESIDUES in line:
                    pdb_data["modelled_residues"] = eval(line.strip(MODELLED_RESIDUES))
                    break
            pdb_data["atoms"] = [ line.strip() for line in F if line.startswith(ATOM)]

            res_list = self._count_residues(pdb_data["atoms"])
            modelled_residues = {"nter" : [], "loops" : [], "cter" : []}

            if len(pdb_data["modelled_residues"]) > 0 and pdb_data["modelled_residues"][0][0] == res_list[0]:
                modelled_residues["nter"].append(pdb_data["modelled_residues"][0])
                pdb_data["modelled_residues"].pop(0)

            if len(pdb_data["modelled_residues"]) > 0 and pdb_data["modelled_residues"][-1][-1] == res_list[-1]:
                modelled_residues["cter"].append(pdb_data["modelled_residues"][-1])
                pdb_data["modelled_residues"].pop(-1)

            modelled_residues["loops"] = pdb_data["modelled_residues"]
            pdb_data["modelled_residues"] = modelled_residues

            pdb_data["residues"] = res_list
        return pdb_data

    def _count_residues(self, atom_list):
        res_list = set([ int(atom[22:30]) for atom in atom_list if atom.startswith(ATOM) ])
        res_list = sorted(res_list)
        return res_list

    def _count_atoms(self, atom_list):
        count = 0
        for atom in atom_list:
            if ATOM in atom:
                count += 1
        return count

    def _trim_residues(self, original_modelled_residues, ter_cutoff, loop_cutoff):
        modelled_residues  = copy.deepcopy(original_modelled_residues)

        for begin, end in modelled_residues["nter"]:
            if end - begin + 1 > ter_cutoff:
                modelled_residues["nter"] = [(end - ter_cutoff + 1, end)]

        for begin, end in modelled_residues["cter"]:
            if end - begin + 1 > ter_cutoff:
                modelled_residues["cter"] = [(begin, begin + ter_cutoff - 1)]

        missing_loops = []
        modelled_loops = []
        for begin, end in modelled_residues["loops"]:
            if end - begin + 1 > loop_cutoff:
                 missing_loops.append((begin, end))
            else:
                modelled_loops.append((begin, end))
        modelled_residues["loops"] = modelled_loops
        modelled_residues["missing_loops"] = missing_loops
        return modelled_residues

    def _trim_atoms(self, modelled_residues, atoms, residue_list):
        """
        :param modelled_residues: dic with keys "nter", "cter", "loops" and "missing_loops"
        :param atoms: list of ATOM lines in pdb
        :param residue_list: list of str
        :return:
        """
        missing_res = []
        for missing in modelled_residues["missing_loops"]:
            missing_res.extend(range(missing[0], missing[1]+1))

        if len(modelled_residues["nter"]) == 1:
            first_res_id = modelled_residues["nter"][0][0]
        else:
            first_res_id = 1

        if len(modelled_residues["cter"]) == 1:
            last_res_id  = modelled_residues["cter"][0][1]
        else:
            last_res_id  = len(residue_list)

        trimed_atoms = []
        for line in atoms:
            resid = int(line[22:30])
            if (first_res_id <= resid <= last_res_id) and (resid not in missing_res):
                trimed_atoms.append(line)

        trimed_atoms = self._insert_ter(trimed_atoms)
        return trimed_atoms

    def _insert_ter(self, atoms):
        """
        insert a "TER" if resid not continuous
        """
        new_atoms = []
        for i in range( len(atoms) - 1):
            new_atoms.append(atoms[i])

            current_resid = int(atoms[i][22:30])
            next_resid    = int(atoms[i+1][22:30])
            if (current_resid != next_resid) and (current_resid+1 < next_resid):
                new_atoms.append("TER")

        new_atoms.append(atoms[-1])
        return new_atoms

    def _residues_to_minimize(self, modelled_residues):
        keys = [key for key in modelled_residues.keys() if key != "missing_loops"]
        residues = []
        for key in keys:
            for begin, end in modelled_residues[key]:
                residues.extend( range(begin, end+1) )
        residues = sorted(residues)
        if len(modelled_residues["nter"]) == 1:
            if modelled_residues["nter"][0][0] != 1:
                shift = 1 - modelled_residues["nter"][0][0]
                residues = [r + shift for r in residues]
        return residues

    def _change_resid(self, atom, new_resid):
        if ATOM not in atom:
            return atom
        entries = atom.split(atom[22:30])
        new_atom = entries[0] + "%4d"%new_resid + " "*4 + entries[1]
        return new_atom

    def _change_resid_sequentially(self):
        atoms = []
        resid = 1
        nlines = len(self._combined_pdb_data["atoms"])
        for i in range( nlines ):
            this_line = self._combined_pdb_data["atoms"][i]
            if ATOM not in this_line:
                atoms.append(this_line)
            else:
                atoms.append( self._change_resid(this_line, resid) )
                this_resid = int(this_line[22:30])
                if i < nlines-1:           # not the last
                    next_line = self._combined_pdb_data["atoms"][i+1]
                    if ATOM in next_line:
                        next_resid = int(next_line[22:30])
                        if next_resid > this_resid:
                            resid += 1
                    else:
                        resid += 1
                else:
                    atoms.append(self._change_resid(this_line, resid))
        return atoms
    
    def _load_ions_cofactors(self, ions_cofactors_files):
        """
        :param ions_cofactors_files: list of str
        :return: None if ions_cofactors_files is empty
        """
        if len(ions_cofactors_files) == 0:
            return None

        for file in ions_cofactors_files:
            pdb_id = os.path.basename(file)[:4]
            if pdb_id != self._pdb_id:
                raise RuntimeError("%s is not from the same pdb with id %s"%(file, self._pdb_id))

        chains = [os.path.basename(file)[4] for file in ions_cofactors_files]
        if len(set(chains).intersection(self._chains)) == 0:
            return None

        pdb_data = {}
        for file in ions_cofactors_files:
            chain = os.path.basename(file)[4]
            if chain in self._chains:
                pdb_data[chain] = self._load_ions_cofactors_file(file)
        return pdb_data

    def _load_ions_cofactors_file(self, file):
        pdb_data = {}
        with open(file, "r") as F:
            pdb_data["atoms"] = [line.strip() for line in F if (line.startswith(HETATM) or line.startswith(TER))]
        pdb_data["natoms"]   = len([line for line in pdb_data["atoms"] if line.startswith(HETATM)])
        pdb_data["nresidues"] = len([line for line in pdb_data["atoms"] if line.startswith(TER)])
        return pdb_data


def parse_modelling_dir(complex_id, modeller_dir):
    """
    :param complex_id: tuple of (pdb_id, chains1, chains2)
    :param modeller_dir: str
    :return:
    """
    pdb_id, chains1, chains2 = complex_id
    modelling_dir = os.path.join(modeller_dir, pdb_id)
    if not os.path.isdir(modelling_dir):
        raise RuntimeError("%s does not exist"%modelling_dir)
    modelled_pdbs = glob.glob( os.path.join(modelling_dir, "*_modelled.pdb") )
    modelled_pdbs = [os.path.basename(file) for file in modelled_pdbs]

    if pdb_id != modelled_pdbs[0][:4]:
        raise RuntimeError("%s does not exist"%pdb_id)
    all_chains = set([pdb[4] for pdb in modelled_pdbs])

    mod_chains1 = False
    for c in chains1:
        if c not in all_chains:
            print("chain %s does not exist in %s"%(c, pdb_id))
            mod_chains1 = True
    if mod_chains1:
        chains1 = [c for c in chains1 if c in all_chains]

    mod_chains2 = False
    for c in chains2:
        if c not in all_chains:
            print("chain %s does not exist in %s"%(c, pdb_id))
            mod_chains2 = True
    if mod_chains2:
        chains2 = [c for c in chains2 if c in all_chains]

    return pdb_id, chains1, chains2, modelling_dir


LIGAND_OUT = "ligand_modelled.pdb"
RECEPTOR_OUT = "receptor_modelled.pdb"

LIGAND_RES_MINIMIZE = "ligand_minimize_list.dat"
RECEPTOR_RES_MINIMIZE = "receptor_minimize_list.dat"


def write_b_receptor_ligand_pdbs(complexes, modeller_dir, ions_cofactors_dir, ter_cutoff=10, loop_cutoff=20):
    """
    complexes:  dict returned by _affinity_data.AffinityData.get_bound_complexes
    """
    print("Combinning chains to form ligands and receptors for ...")
    for name, complex_id in complexes.items():
        print(name)
        if not os.path.isdir(name):
            os.makedirs(name)
        pdb_id, chains1, chains2, modelling_dir = parse_modelling_dir(complex_id, modeller_dir)

        ic_dir = os.path.join(ions_cofactors_dir, name)
        if not os.path.isdir(ic_dir):
            ions_cofactors_files = []
        else:
            ions_cofactors_files = glob.glob(os.path.join(ic_dir, pdb_id + "*" + ".pdb"))

        partners = [ ChainCombine(pdb_id, chains, modelling_dir, ions_cofactors_files) for chains in (chains1, chains2) ]
        for p in partners:
            p.trim_residues(ter_cutoff, loop_cutoff)
            p.combine()
        partners.sort(key=lambda c: c.get_nresidues())
        partners[0].write_pdb(out = os.path.join(name, LIGAND_OUT))
        partners[0].write_residues_to_minimize(os.path.join(name, LIGAND_RES_MINIMIZE))

        partners[1].write_pdb(out = os.path.join(name, RECEPTOR_OUT))
        partners[1].write_residues_to_minimize(os.path.join(name, RECEPTOR_RES_MINIMIZE))

    print("Done combinning chains")
    print("")
    return None

