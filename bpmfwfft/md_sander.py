"""
python wraper of sander (AmberTools) to calculate energy 
"""
from __future__ import print_function

import os

import numpy as np

SANDER_PHASES = ["sander_Gas", "sander_PBSA", "sander_HCT", "sander_OBC1", "sander_OBC2", "sander_GBn", "sander_GBn2"] 


def _sander_script_gen(phase, no_bonded):
    """
    adapted from AlGDock _sander_energy
    :param phase: str, one of SANDER_PHASES
    :param no_bonded: bool
    :return: str, which can be written to file and pass to "-i" option of sander
    """
    assert phase in SANDER_PHASES, "Phase %s is unknown"%phase 

    script = """Calculating energies with sander
&cntrl
    imin=5,    ! read trajectory in for analysis
    maxcyc=1,  ! Single-point energy calculation on each frame
    ntx=1,     ! input is read formatted with no velocities
    irest=0,   ! do not restart
    ntb=0,     ! no periodicity and no PME
    idecomp=0, ! no decomposition
    ntc=1,     ! No SHAKE
    cut=9999., !
    """
    if no_bonded:
        script += """ntf=7,     ! No bond, angle, dihedral forces calculated"""
    else:
        script += """ntf=1,     ! complete interaction calculation"""

    if phase=='sander_Gas':
        script += """
/
        """

    elif phase=='sander_PBSA':
        script += """
    ipb=2,     ! Default PB dielectric model
    inp=2,     ! non-polar from cavity + dispersion
/
&pb
    radiopt=0, ! Use atomic radii from the prmtop file
    fillratio=4.0,
    sprob=1.4,
    cavity_surften=0.0378, ! (kcal/mol) Default in MMPBSA.py
    cavity_offset=-0.5692, ! (kcal/mol) Default in MMPBSA.py
/
        """
    else:
        key = phase.split('_')[-1]
        igb = {'HCT':1, 'OBC1':2, 'OBC2':5, 'GBn':7, 'GBn2':8}[key]
        script += """
    igb=%d,     !
    gbsa=2,    ! recursive surface area algorithm (for postprocessing)
/
        """%(igb)
    return script


def _array_2_mdcrd(crd):
    """
    :param crd: ndarray of shape (nconfs, natoms, 3)
    :return: str, which can be written to mdcrd file
    """
    assert len(crd.shape) == 3, "crd must be of 3D"
    assert crd.shape[2] == 3, "crd last dimension must be 3"  

    out_str = "\n"
    for conf in crd:
        count = 0
        for atom in conf:
            for xyz in atom:
                count += 1
                out_str += "%8.3f"%(xyz)
                if (count%10) == 0:
                    out_str += "\n"
        out_str += "\n"
    return out_str


def _array_2_inpcrd(crd):
    """
    :param crd: ndarray with shape (natoms, 3)
    :return: str, which can be written to inpcrd file
    """
    assert len(crd.shape) == 2, "crd must be of 2D"
    assert crd.shape[1] == 3, "crd last dimension must be 3"

    out_str = """default_name\n"""
    out_str += """%6d\n"""%crd.shape[0]
    for i in range(crd.shape[0]):
        x, y, z = crd[i]
        out_str += """%12.7f%12.7f%12.7f"""%(x, y, z)
        if (i+1)%2 == 0:
            out_str += """\n"""
    return out_str


def _extract_component_energies(file):
    exclude_terms = ["VDWAALS", "1-4VDW", "EEL", "1-4EEL"]
    with open(file, "r") as handle:
        data = handle.read().strip().split(" BOND")
    data.pop(0)
    assert len(data) > 0, file + " has no energy output"
    energies = []
    for rec in data:
        entries = rec[:rec.find("\nminimization")].replace('1-4 ','1-4').split()
        terms = entries[2::3]
        e = 0.
        for i in range(len(entries)):
            if (entries[i] in terms) and (entries[i] not in exclude_terms):
                e += np.float(entries[i+2])
        energies.append(e)
    return np.array(energies, dtype=float)

def _extract_total_energies(file):
    energies = []
    out_file_lines = open(file, "r").readlines()
    for i in range(len(out_file_lines)):
        line = out_file_lines[i]
        if ("ENERGY" in line) and (line.strip().startswith("NSTEP")):
            e = out_file_lines[i+1].split()[1]
            e = np.float(e)
            energies.append(e)
    return np.array(energies, dtype=float)

def sander_energy(prmtop_file, crd, phase, tmp_dir, no_bonded=True):
    """
    crd in angstroms
    crd:    np.ndarray with shape  (nconfs, natom, 3)
    """
    if len(crd.shape) == 2:
        crd = np.array([crd], dtype=float)
    nconfs = crd.shape[0]

    tmp_dir = os.path.abspath(tmp_dir)
    sander_file_prefix = "sander_tmp_%d"%np.random.randint(0,10000)

    SANDER_SCRIPT = os.path.join(tmp_dir, sander_file_prefix+".inp")
    INPCRD = os.path.join(tmp_dir, sander_file_prefix+".inpcrd")
    SANDER_OUT = os.path.join(tmp_dir, sander_file_prefix+".out")
    MDCRD_FILE = os.path.join(tmp_dir, sander_file_prefix+".mdcrd")
    RESTART = os.path.join(tmp_dir, sander_file_prefix+".restrt")

    script_str = _sander_script_gen(phase, no_bonded)
    open(SANDER_SCRIPT, "w").write(script_str)

    crd_str = _array_2_mdcrd(crd)
    open(MDCRD_FILE, "w").write(crd_str)
    del crd_str

    inpcrd_str = _array_2_inpcrd(crd[0])
    open(INPCRD, "w").write(inpcrd_str)
    del inpcrd_str
    del crd

    #run_sanser = "sander -O -i " + SANDER_SCRIPT + " -o " + SANDER_OUT + " -p " + prmtop_file + \
    #        " -c " + INPCRD + " -y " + MDCRD_FILE + " -r " + RESTART
    #print run_sanser
    #os.system(run_sanser)

    import subprocess
    args_list = ["sander", "-O", "-i", SANDER_SCRIPT, "-o", SANDER_OUT, "-p", prmtop_file, "-c", INPCRD, "-y", MDCRD_FILE, "-r", RESTART]
    print(' '.join(args_list))
    p = subprocess.Popen(args_list)
    p.wait()

    #remove_tmp_files = "rm " + " ".join([SANDER_SCRIPT, SANDER_OUT, INPCRD, MDCRD_FILE, RESTART])
    #remove_tmp_files = "rm " + " ".join([SANDER_OUT, RESTART])
    #print remove_tmp_files
    #os.system(remove_tmp_files)
    if no_bonded:
        energies = _extract_component_energies(SANDER_OUT)
    else:
        energies = _extract_total_energies(SANDER_OUT)

    if len(energies) != nconfs:
        raise RuntimeError("The number of energy values is not the same as number of confs")
    return energies


