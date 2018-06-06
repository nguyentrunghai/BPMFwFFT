
"""
define class to handle the affinity data pushlished at
http://bmm.crick.ac.uk/~bmmadmin/Affinity/
"""

import csv
import numpy as np
import pandas as pd


class AffinityData(object):
    """
    load the affinity data files
    give interface to access the data
    """
    def __init__(self, affinity_data_files):
        """
        :param affinity_data_files: list of str
        """
        self._data_frame = self._load_tsvfiles(affinity_data_files)

    def unique_pdb_ids(self):
        """
        return a set of unique pdb ids to be downloaded
        """
        ids = []
        for col in ["Complex PDB", "Unbound PDB Protein A", "Unbound PDB Protein B"]:
            ids.extend( list(self._data_frame[col].values) )
        ids = [id.split("_")[0].lower() for id in ids]
        return set(ids)

    def get_bound_complexes(self):
        """
        Each name is unique, and can be used to identify the complex.
        The structure of names is "XXXX_AB:FGH", where "XXXX" is 4-letter pdb id,
        "AB" is chains of protein A, "FGH" is chains of protein B.

        return a dic { name:(pdb_id (str), chains1, chains2 ) }
                    name: str; chains1:  list of str; chains2:  list of str
        """
        names = list(self._data_frame["Complex PDB"])
        complex_chains = {}
        for name in names:
            pdb_id, chains = name.split("_")
            pdb_id = pdb_id.lower()
            chains1, chains2 = chains.split(":")
            chains1 = [c for c in chains1]
            chains2 = [c for c in chains2]
            if len(chains1) == 0 or len(chains2) == 0:
                raise RuntimeError("%s does not have one or both binding partners"%name)
            complex_chains[name] = (pdb_id, chains1, chains2)
        return complex_chains

    def get_col_names(self):
        return list(self._data_frame.columns)

    def get_data_from_col(self, col_name):
        complex_names = list(self._data_frame["Complex PDB"])
        col = {}
        for name in complex_names:
            row = self._data_frame["Complex PDB"] == name
            value = self._data_frame[col_name][row].values
            if value.shape != (1,):
                raise RuntimeError("There are more than one %s at %s"%(col_name, name))
            col[name] = value[0]
        return col

    def get_dG(self):
        """
        return a dic, dG[complex_name] -> dG
        """
        dG = self.get_data_from_col("dG")
        for name in dG.keys():
            dG[name] = np.float(dG[name])
        return dG

    def _load_tsvfile(self, file):
        """
        load tsv file
        file is a str
        return a pandas.DataFrame object
        TODO: use pd.read_table
        """
        with open(file) as tsvfile:
            tsvreader = csv.reader(tsvfile, delimiter="\t")
            tsvreader.next()                            # ignore the first line
            data_fields = tsvreader.next()
            # in the the file, there are two columns with the same name "Unbound PDB"
            # make each field in data_fields unique
            for i in range(len(data_fields)):
                if data_fields[i] == "Unbound PDB":
                    data_fields[i] += " " + data_fields[i+1]

            # records is a list of dict
            records = []
            for line in tsvreader:
                if len(line) != len(data_fields):
                    raise RuntimeError("line %s does not have the same len as %s" %("\t".join(line), "\t".join(data_fields)))
                # put each line into a dict
                tmp = {}
                for i in range( len(data_fields) ):
                    tmp[data_fields[i]] = line[i]
                records.append(tmp)
        return pd.DataFrame(records)
    
    def _load_tsvfiles(self, files):
        """
        load multiple tsv files
        return a concatenated pd.DataFrame
        files is a list of str
        """
        assert len(files) == len(set(files)), "some files have the same name"
        frames = [self._load_tsvfile(file) for file in files]
        return pd.concat(frames)
    


