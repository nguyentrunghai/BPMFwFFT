"""
to extract bound ions and cofactors from pdb files
"""
from __future__ import print_function

import os
import argparse

from _affinity_data import AffinityData

parser = argparse.ArgumentParser()
parser.add_argument("--affinity_dir", type=str, default="affinity")

# this contains files downloaded from pdb, not modelled by modeller
parser.add_argument("--pdb_dir", type=str, default="download_pdb")

parser.add_argument("--include_all", action="store_true", default=False)

args = parser.parse_args()

AFFINITY_DATA_FILES = ["affinity_v1.tsv",  "affinity_v2.tsv"]
AFFINITY_DATA_FILES = [os.path.join(args.affinity_dir, file) for file in AFFINITY_DATA_FILES]

CHAIN_POS = (21, 22)
RES_NAME_POS = (17, 20)
RES_ID_POS = (22, 26)
ATOM_NAME_POS = (12, 17)

ACCEPTED_IONS = ["MG", "CA", "FE", "ZN", "SR"]
ACCEPTED_COFACTORS = ["ADP", "ATP", "GDP", "GTP"]


def _extract_ions_cofactors_from_pdb(pdb_file, chains, ions_cofactors, out_dir):
    """
    pdb_file:   str, name of original pdb file
    chains:     list of str, list of chains for which ions and/or cofactors are extracted
    ions_cofactors: list of str, list of cofactors name

    for each chain, the function will write coordinates of ions and cofactors to 
    id+chain+"_ions_cofactors.pdb"
    """
    id = os.path.basename(pdb_file)[:-4]
    all_hetatm = [line for line in open(pdb_file, "r") if line.startswith("HETATM")]
    if len(all_hetatm) == 0:
        print(pdb_file + " does not have any HETATM")
        return None

    i, j = CHAIN_POS
    all_hetatm = [line for line in all_hetatm if line[i:j] in chains]
    if len(all_hetatm) == 0:
        print(pdb_file + " does not have any HETATM with chain in ", chains)
        return None

    i, j = RES_NAME_POS
    all_hetatm = [line for line in all_hetatm if line[i:j].strip() in ions_cofactors]
    if len(all_hetatm) == 0:
        print(pdb_file + " does not have any HETATM with resname in ", ions_cofactors)
        return None

    atom_records = {}
    for i, line in enumerate(all_hetatm):
        current_resid = int(line[RES_ID_POS[0] : RES_ID_POS[1]])

        chain = line[CHAIN_POS[0] : CHAIN_POS[1]]
        if chain not in atom_records.keys():
            atom_records[chain] = {}

        resname = line[RES_NAME_POS[0] : RES_NAME_POS[1]].strip()
        try:
            len(atom_records[chain][resname])
        except KeyError:
            atom_records[chain][resname] = []

        atom_records[chain][resname].append(line)

        if i+1 < len(all_hetatm):
            next_resid = int(all_hetatm[i+1][RES_ID_POS[0] : RES_ID_POS[1]])
        else:
            next_resid = current_resid
            atom_records[chain][resname].append("TER\n")

        if current_resid != next_resid:
            atom_records[chain][resname].append("TER\n")

    if not os.path.isdir(out_dir):
        os.system("mkdir "+out_dir)

    for chain in atom_records.keys():
        with open(os.path.join(out_dir, id + chain + ".pdb"), "w") as F:
            for resname in atom_records[chain].keys():
                F.write("".join([line for line in atom_records[chain][resname]]))
    return None


def extract_ions_cofactors_from_database(acepted_ions, acepted_cofactors, 
                                        original_pdb_dir, include_all):
    """
    acepted_resnames:   list of str
    return dict  sel[complex_name] -> list of ions and cofactors' name
    """
    aff = AffinityData(AFFINITY_DATA_FILES)
    ions_cofactors_bound = aff.get_data_from_col("Cofactors (bound)")
    ions_cofactors_binding = aff.get_data_from_col("Cofactors (binding)")

    acepted_residues = acepted_ions + acepted_cofactors

    sel = {}
    if include_all:
        for complex in ions_cofactors_bound.keys():
            pdb_id = complex[:4].lower()
            chains = complex.split("_")[-1]
            chains = [c for c in chains if c != ":"]
            with open(os.path.join(original_pdb_dir, pdb_id+".pdb"), "r") as F:
                for line in F:
                    if line.startswith("HETATM"):
                        chain = line[CHAIN_POS[0] : CHAIN_POS[1]]
                        res = line[RES_NAME_POS[0] : RES_NAME_POS[1]].strip()
                        if (chain in chains) and (res in acepted_residues):
                            if complex not in sel.keys():
                                sel[complex] = []
                            sel[complex].append(res)

    for complex in ions_cofactors_bound.keys():
        a = ions_cofactors_bound[complex].split(",")
        for i in a:
            b = i.split("/")
            for j in b:
                if j.strip().upper() in acepted_residues:
                    if complex not in sel.keys():
                        sel[complex] = []
                    sel[complex].append(j.upper())

    for complex in ions_cofactors_binding.keys():
        a = ions_cofactors_binding[complex].split(",")
        for i in a:
            b = i.split("/")
            for j in b:
                if j.strip().upper() in acepted_residues:
                    if complex not in sel.keys():
                        sel[complex] = []
                    sel[complex].append(j.upper())

    for complex in sel.keys():
        sel[complex] = list(set(sel[complex]))

    return sel


ions_cofactors = extract_ions_cofactors_from_database(ACCEPTED_IONS, ACCEPTED_COFACTORS,
                                                        args.pdb_dir, include_all=args.include_all)
for complex in ions_cofactors:
    print(complex, ions_cofactors[complex])
    pdb_id = complex[:4].lower()
    chains = complex.split("_")[-1]
    chains = [c for c in chains if c != ":"]
    print(chains)

    original_pdb = os.path.join(args.pdb_dir, pdb_id+".pdb")
    _extract_ions_cofactors_from_pdb(original_pdb, chains, ions_cofactors[complex], complex)

print("DONE")
