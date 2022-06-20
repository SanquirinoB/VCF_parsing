from itertools import product
from Structures import MetaRef
import os
import shutil

class ReferenceProcessor():
    def __init__(self, _raw_reference_path, _destination_folder):
        # First parameters
        self.n_refs = 0
        self.refs_len = 0
        self.raw_reference_path = _raw_reference_path

        # We asume that VCFParser is excecuted first
        tmp_folder = os.path.join(_destination_folder, "Tmp")
        self.parsing_folder = os.path.join(tmp_folder, "Parsing")
        self.metadata_folder = os.path.join(tmp_folder, "Meta_data")

        # Helper variables
        self.meta_structure = MetaRef()
        self.current_ref_data = {}
        self.reference_data = {}
        self.checkpoint_refs_len = 0

        self.raw_ref_file = None
        self.meta_file = None

    def IsDescriptionLine(self, line):
        return line[0] == ">"

    def GenerateCharacterization(self, line):
        values = line[1:].split(" ")
        # Is mandatory that the first element will be the id
        self.current_ref_data["ID"] = values[0]
        self.current_ref_data["internal_ID"] = self.n_refs
        self.current_ref_data["n_bases"] = 0
        self.n_refs += 1

    def SaveRefData(self):
        self.meta_structure.m_ID = self.current_ref_data["internal_ID"]
        self.meta_structure.m_nBases = self.current_ref_data["n_bases"]
        self.meta_structure.m_relPos = self.checkpoint_refs_len
        self.meta_file.write(bytearray(self.meta_structure))

        # For parsing use
        # ('Y', {'length': '59373566', 'assembly': 'b37', 'relPosRef': 3042602996, 'internalID': 84, 'ID': 'Y'})
        dict_aux = {}
        dict_aux["ID"] = self.current_ref_data["ID"]
        dict_aux["internal_ID"] = self.current_ref_data["internal_ID"]
        dict_aux["length"] = self.current_ref_data["n_bases"]
        dict_aux["relPosRef"] = self.checkpoint_refs_len

        self.reference_data[self.current_ref_data["ID"]] = dict_aux.copy()       

    def AppendBasePermutations(self, processed_reference):
        bases = ["A", "C", "T", "G", "N"]
        edits = []
        for i in range(1,5):
            for p in product(bases, repeat=i):
                edits.append("".join(p))

        for edit in edits:
            dict_aux = {}
            dict_aux["ID"] = edit
            dict_aux["internal_ID"] = self.n_refs
            dict_aux["length"] = len(edit)
            dict_aux["relPosRef"] = self.refs_len
            processed_reference.write(edit)

            self.n_refs += 1
            self.refs_len += len(edit)

    def StartReferenceProcessing(self):
        assert os.path.isdir(self.parsing_folder)
        assert os.path.isdir(self.metadata_folder)

        # Get propper paths
        processed_ref_file_path = os.path.join(self.parsing_folder, "Reference.tmprlz")
        metadata_file_path = os.path.join(self.metadata_folder, "Reference.metarlz")

        # Open metadata file for refs info
        self.meta_file = open(metadata_file_path, 'wb')

        # Before we process the first ref, we start to check if all the values are correct
        start_checking = False
        self.checkpoint_refs_len = 0

        # Start processing
        with open(self.raw_reference_path, 'r') as raw_reference, open(processed_ref_file_path, 'w') as processed_reference:
            # First we read the first line of the ref, which should be like >1 dna:chromosome chromosome:GRCh37:1:1:249250621:1\n
            line = raw_reference.readline().rstrip()
            while line:
                if(self.IsDescriptionLine(line)):
                    # Before we start to check a new ref, we need to ensure all info read before was correct
                    if(start_checking):
                        assert (self.refs_len - self.checkpoint_refs_len) == self.current_ref_data["n_bases"]
                        self.SaveRefData()

                    self.checkpoint_refs_len = self.refs_len
                    self.GenerateCharacterization(line)
                    start_checking = True
                else:
                    self.current_ref_data["n_bases"] += len(line)
                    self.refs_len += len(line)
                    processed_reference.write(line)
                line = raw_reference.readline().rstrip()

            self.AppendBasePermutations(processed_reference)
            

        self.meta_file.close()

    def GetReferenceData(self):
        return self.refs_len, self.reference_data


            
        

