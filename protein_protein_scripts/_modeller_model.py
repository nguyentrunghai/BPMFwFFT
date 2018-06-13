"""
"""
import glob

from modeller import environ
from modeller import selection
from modeller.automodel import automodel
from modeller.automodel import loopmodel
from modeller.automodel import refine


class ModellingError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)


def run_modeller(id, chain, missing_res_ranges, sequence_len, align_file):
    """
    :param id:
    :param chain:
    :param missing_res_ranges:
    :param sequence_len:
    :param align_file:
    :return:
    """
    env = environ()
    env.io.atom_files_directory = ['.']
    
    modifying_string = _redefine_loopmodel_string(missing_res_ranges, sequence_len)

    exec(modifying_string)
    model = MyModel(env, alnfile=align_file, knowns=(id), sequence = id+chain+"_full")
    model.make()

    _copy_result(id, chain, missing_res_ranges)
    return None


def _redefine_loopmodel_string(missing_res_ranges, sequence_len):
    """
    retrun a str to be run by exec to modyfy the class loopmodel.
    :param missing_res_ranges:
    :param sequence_len:
    :return:
    """
    out_string = """class MyModel(automodel):\n\tdef select_atoms(self):\n\t\treturn selection( """

    if len(missing_res_ranges) == 0:
        out_string += "self.residue_range(0, 0) )\n"
    else:
        for i in range(len(missing_res_ranges)-1):
            out_string += "self.residue_range(%d, %d), "%missing_res_ranges[i]
        out_string += "self.residue_range(%d, %d) )\n"%missing_res_ranges[-1]

    return out_string


def _copy_result(id, chain, missing_res_ranges):
    """
    :param id:
    :param chain:
    :param missing_res_ranges:
    :return:
    """
    missing_res_ranges = [(star+1, end+1) for star, end in missing_res_ranges]
    output_pdbs = glob.glob(id+chain+"_full.B*.pdb")
    if len(output_pdbs) > 0:
        pdb_text = open(output_pdbs[0], "r").read()
        remark_str = "REMARK  MODELLED RESIDUES: " + str(missing_res_ranges) + "\n"
        pdb_text = remark_str + pdb_text

        open(id+chain+"_modelled.pdb", "w").write(pdb_text)
    return None


