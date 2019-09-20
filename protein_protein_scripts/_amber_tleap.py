"""
functions to generate AMBER topology and coordinates files
"""

from __future__ import print_function

import os
import glob


def _parse_cofactors_prep_dir(dir):
    """
    return a dic [cofactor name "ATP"] -> (prep_file, frcmod_file)
    """
    prep_files = glob.glob(os.path.join(dir, "*.prep"))
    prep_files = [os.path.abspath(file) for file in prep_files]

    frcmod_files = glob.glob(os.path.join(dir, "*.frcmod")) 
    frcmod_files = [os.path.abspath(file) for file in frcmod_files]

    cofactor_names = [os.path.basename(file)[:-5] for file in prep_files]
    print(cofactor_names)

    cofactors_prep = {}
    for name in cofactor_names:
        for file in prep_files:
            if os.path.basename(file).startswith(name):
                prep_file = file
        
        for file in frcmod_files:
            if os.path.basename(file).startswith(name):
                frcmod_file = file
        cofactors_prep[name] = (prep_file, frcmod_file)
    return cofactors_prep


def _parse_cofactors_in_pdb(pdb_file, allowed_cofactor_names):
    """
    pdb_file:   str
    allowed_cofactor_names: list of str
    """
    RES_NAME_POS  = (17, 20)
    cofactor_names = []
    with open(pdb_file, "r") as F:
        for line in F:
            if line.startswith("HETATM"):
                resname = line[RES_NAME_POS[0] : RES_NAME_POS[1]].strip()
                if resname in allowed_cofactor_names:
                    cofactor_names.append(resname)
    return list(set(cofactor_names))


def _write_tleap_script(ligand_pdb, receptor_pdb, cofactors_prep, out_dir, 
        ligand_prefix, receptor_prefix, complex_prefix, tleap_script):
    """
    ligand_pdb: str
    receptor_pdb:   str
    cofactors_prep: dic returned by _parse_cofactors_prep_dir 
    """
    out_dir = os.path.abspath(out_dir)

    allowed_cofactor_names = cofactors_prep.keys()
    cofactor_names_lig = _parse_cofactors_in_pdb(os.path.join(out_dir, ligand_pdb), 
                                                    allowed_cofactor_names)
    cofactor_names_rec = _parse_cofactors_in_pdb(os.path.join(out_dir, receptor_pdb), 
                                                    allowed_cofactor_names)
    cofactor_names = set(cofactor_names_lig + cofactor_names_rec)

    script = "\n"
    script += "source leaprc.ff14SB\n"
    script += "set default PBRadii mbondi2\n"
    script += "loadamberparams frcmod.ionslm_1264_tip3p\n"

    if len(cofactor_names) > 0:
        for name in cofactor_names:
            script += "loadamberprep " + cofactors_prep[name][0] + "\n"
            script += "loadamberparams " + cofactors_prep[name][1] + "\n"

    script += "receptor = loadpdb %s\n"%receptor_pdb
    script += "ligand = loadpdb %s\n\n"%ligand_pdb

    lig_prmtop = ligand_prefix + ".prmtop"
    lig_inpcrd = ligand_prefix + ".inpcrd"
    lig_pdb = ligand_prefix + ".pdb"
    script += "saveAmberParm ligand %s %s\n"%(lig_prmtop, lig_inpcrd)
    script += "savepdb ligand %s\n\n"%lig_pdb

    rec_prmtop = receptor_prefix + ".prmtop"
    rec_inpcrd = receptor_prefix + ".inpcrd"
    rec_pdb    = receptor_prefix + ".pdb"
    script += "saveAmberParm receptor %s %s\n"%(rec_prmtop, rec_inpcrd)
    script += "savepdb receptor %s\n\n"%rec_pdb

    script += "complex = combine {receptor ligand}\n\n"

    comp_prmtop = complex_prefix + ".prmtop"
    comp_inpcrd = complex_prefix + ".inpcrd"
    comp_pdb = complex_prefix + ".pdb"
    script += "saveAmberParm complex %s %s\n"%(comp_prmtop, comp_inpcrd)
    script += "savepdb complex %s\n\n"%comp_pdb

    script += "quit\n"

    script_file = os.path.join(out_dir, tleap_script)
    open(script_file, "w").write(script)
    return None


TLEAP_FAILED = "TLEAP_FAILED"



def _run_tleap(tleap_script_file):
    if not os.path.isfile(tleap_script_file):
        raise RuntimeError("%s does not exist" % tleap_script_file)

    tleap_script = os.path.abspath(tleap_script_file)
    cwd = os.getcwd()
    run_dir = os.path.dirname(tleap_script)
    script_name = os.path.basename(tleap_script)

    os.chdir(run_dir)
    print("\n\n\n\nRunning tleap for " + tleap_script_file)
    print("cwd", os.getcwd())
    os.system("tleap -f %s" % script_name)

    # check results

    prmtop_files = glob.glob("*.prmtop")
    if len(prmtop_files) != 3:
        os.system("touch "+TLEAP_FAILED)
    else:
        for file in prmtop_files:
            if not os.path.exists(file):
                os.system("touch "+TLEAP_FAILED)
                break
            else:
                if os.path.getsize(file) == 0:
                    os.system("touch "+TLEAP_FAILED)
                    break
    os.chdir(cwd)
    return None


LIGAND_PDB_INP = "ligand_modelled.pdb"
RECEPTOR_PDB_INP = "receptor_modelled.pdb"

LIGAND_OUT_PREFIX = "ligand"
RECEPTOR_OUT_PREFIX = "receptor"
COMPLEX_OUT_PREFIX = "complex"
TLEAP = "setup.tleap"


def generate_prmtop(cofactors_frcmod_dir):
    """
    :param cofactors_prep_dir: str
    :return: None
    """
    complex_names = glob.glob("*")
    complex_names = [os.path.basename(d) for d in complex_names if os.path.isdir(d)]

    cofactors_prep = _parse_cofactors_prep_dir(cofactors_frcmod_dir)

    if len(complex_names) == 0:
        print("Do nothing!")
        return None
    print("Generating amber top for ...")
    for complex_name in complex_names:
        print(complex_name)
        out_dir = os.path.abspath(complex_name)
        _write_tleap_script(LIGAND_PDB_INP, RECEPTOR_PDB_INP, cofactors_prep, out_dir,
                            LIGAND_OUT_PREFIX, RECEPTOR_OUT_PREFIX, COMPLEX_OUT_PREFIX, TLEAP)

        tleap_script_file = os.path.join(out_dir, TLEAP)
        _run_tleap(tleap_script_file)

    print("Done with Amber")
    print("Checking for failed ...")
    for complex_name in complex_names:
        if os.path.exists(os.path.join(complex_name, TLEAP_FAILED)):
            print(complex_name, TLEAP_FAILED)
    return None

