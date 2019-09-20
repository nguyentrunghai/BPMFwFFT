
"""
define class to handle pdb files
"""

import modeller
import modeller.automodel
import Bio.SeqIO

from _modeller_model import run_modeller


class AddMissing(object):
    def __init__(self, pdb_file):
        """
        :param pdb_file: str
        """
        self._file = pdb_file
        self._check_pdb()
        self._id = pdb_file[:-4]
        self._text_lines = open(pdb_file, "r").readlines()

        self._structureX_seq_header = self._structureX_seq_from_modeller()
        self._full_sequences = self._full_seq_from_Bio()
        self._chains_list = self._get_list_of_chains()
        self._residue_ranges = self._search_res_ranges()

    def write_dummy_alignments(self):
        """
        write dummy alignment files for all chains
        """
        for chain, full_seq in self._full_sequences:
            res_start, res_end = self._residue_ranges[chain]
            header = self._modify_structureX_seq_header(res_start, chain, res_end, chain)
            open(self._id + chain +".seq", "w").write(full_seq + header)
        return None

    def do_auto_align(self):
        """
        modeller auto align reads all the dummy alignment files
        to align the full sequence for each chain with the sequence appearing in ATOM
        set self._alignment_files and self._missing_res_ranges
        the indices in self._missing_res_ranges is pythonic but inclusive
        """
        self._alignment_files = {}
        self._missing_res_ranges = {}
        self._seq_lengths = {}

        for chain in self._chains_list:
            dummy_align_file = self._id + chain +".seq"
            align_file  = dummy_align_file + ".ali"

            self._modeller_auto_alignment(chain, dummy_align_file)
            self._check_align_file(chain, align_file)

            self._alignment_files[chain] = align_file
            self._missing_res_ranges[chain] = self._missing_residues(align_file) 
            self._seq_lengths[chain] = self._len_sequence(chain, align_file)
        return None

    def do_auto_modeller(self):
        for chain in self._chains_list:
            run_modeller(self._id, chain, self._missing_res_ranges[chain], 
                    self._seq_lengths[chain], self._alignment_files[chain] )
        return None

    def write_vmd_script(self, file="vmd.tcl"):
        with open(file, "w") as F:
            F.write("mol new " + self._id+".pdb" + "\n")
            for chain in self._chains_list:
                F.write("mol new " + self._id+chain+"_modelled"+".pdb" + "\n")
        return None

    def _check_pdb(self):
        with open(self._file, "r") as pdb:
            pdb_text = pdb.read()
            for rec in ["SEQRES", "ATOM"]:
                if rec not in pdb_text:
                    raise RuntimeError("%s not in %s"%(rec, self._file))
        return None

    def _structureX_seq_from_modeller(self):
        """
        return a str containing the first two lines of the sequence corresponding to structureX 
        a file named [self._id]_structureX.seq also written
        """
        env = modeller.environ()
        model = modeller.model(env, file=self._id)
        aln = modeller.alignment(env)
        aln.append_model(model, align_codes=self._id)
        out_file = self._id+"_structureX.seq"
        aln.write(file = out_file)
        out_str = open(out_file, "r").read()
        out_str = [ c for c in out_str.split("\n") if c]
        out_str = "\n".join(out_str[:2]) + "\n*"
        return out_str

    def _modify_structureX_seq_header(self, resid_start, chain_start, resid_end, chain_end):
        """
        use self._structureX_seq_header as a template to modify
        the header for a given chain
        return a new str
        TODO, force input information to be correct w.r.t. the pdb
        """
        assert type(resid_start) == str and type(chain_start) == str and \
                type(resid_end) == str and type(chain_end) == str, "Input type not correct"
        header_lines = self._structureX_seq_header.split("\n")

        description = header_lines[1].split(":")
        description[2] = " %s "%resid_start
        description[3] = " %s "%chain_start
        description[4] = " %s "%resid_end
        description[5] = " %s "%chain_end
        description = ":".join(description)

        header_lines[1] = description
        return "\n".join(header_lines)

    def _full_seq_from_Bio(self):
        """
        return a list of tuples, each is a (chain : pir_sequence) pair.
        chains are ordered as they appear in the pdb
        """
        records = Bio.SeqIO.parse(self._file, "pdb-seqres")
        records = [rec.format("fasta") for rec in records]

        sequences = []
        for seq in records:
            lines = seq.split("\n")
            chain   = lines[0].split(":")[1][0]
            seq_str  = ">P1;" + self._id + chain + "_full" + "\n"
            seq_str += "sequence:::::::::" + "\n"

            aa_seq_str = "".join(lines[1:])
            if len(aa_seq_str) > 0:
                count = 0
                for c in aa_seq_str:
                    seq_str += c
                    count += 1
                    if count % 75 == 0:
                        seq_str += "\n"
                seq_str += "*\n"
                sequences.append( (chain, seq_str) )
        return sequences

    def _get_list_of_chains(self):
        chains = [seq[0] for seq in self._full_sequences]
        return chains

    def _search_res_ranges(self):
        """
        return a dic with { chain:(res_start, res_end) }
        """
        res_ranges = {}
        for chain in self._chains_list:

            chain_text_lines = []
            for line in self._text_lines:
                if line[0:4] == "ATOM" and line[21:22] == chain:
                    chain_text_lines.append(line)

            if len(chain_text_lines) == 0:
                raise RuntimeError("%s does not have chain %s"%(self._file, chain))

            start = chain_text_lines[0][22:27]
            end   = chain_text_lines[-1][22:27]
            res_ranges[chain] = (start, end)
        return res_ranges

    def _modeller_auto_alignment(self, chain, dummy_align_file):
        """
        """
        env = modeller.environ()
        model = modeller.automodel.automodel(env, alnfile=dummy_align_file, 
                knowns=self._id, sequence = self._id+chain+"_full") 
        model.auto_align()
        return None

    def _load_align_file(self, file):
        """
        return a dict of list, whose keys are sequence id
        """
        seqs = open(file, "r").read()
        seqs = seqs.split(">P1;")
        seqs = [s for s in seqs if s and s != "\n"]

        out_seqs = {}
        for seq in seqs:
            lines = seq.split("\n")
            out_seqs[lines[0]] = [ lines[1], "".join(lines[2:]) ]
        return out_seqs

    def _len_sequence(self, chain, align_file):
        align_seqs = self._load_align_file(align_file)
        if len(align_seqs[self._id + chain + "_full"][-1]) != len(align_seqs[self._id][-1]):
            raise RuntimeError("len of the two seqs in %s not the same"%align_file)
        return len(align_seqs[self._id + chain + "_full"][-1])

    def _check_align_file(self, chain, file):
        """
        raise exception 
                        if [id]_full contains missing "-"
                        if any of aligned sequences has chain break marker "/" 
        """
        align_seqs = self._load_align_file(file)
        if "-" in align_seqs[self._id + chain + "_full"][-1]:
            raise RuntimeError("%s has missing aa in %s"%(file, self._id + chain + "_full"))
        
        for id in align_seqs.keys():
            if "/" in align_seqs[id][-1]:
                raise RuntimeError("%s has '/' in %s"%(file, id))
        return None

    def _missing_residues(self, align_file):
        """
        look for "-" in the align_file
        return list of pairs of ( start, end ), pythonic indexing but inclusive
        """
        seq = self._load_align_file(align_file)[self._id][-1]
        seq = "".join([ c for c in seq if c != "*" ])

        missing_ind = []
        for i in range(len(seq)):
            if seq[i] == "-":
                missing_ind.append(i)

        if len(missing_ind) == 0:
            return []

        missing_res_ranges = []
        begin = missing_ind[0]
        for i in range(len(missing_ind)-1):
            if missing_ind[i] != missing_ind[i+1]-1:
                end = missing_ind[i]
                missing_res_ranges.append((begin, end))
                begin = missing_ind[i+1]
        end = missing_ind[-1]
        missing_res_ranges.append((begin, end))
        return missing_res_ranges


