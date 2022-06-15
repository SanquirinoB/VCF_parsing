from itertools import permutations
from operator import truediv
import os
import shutil

class ReferenceProcessor():
    def __init__(self, _n_refs, _refs_len, _bases_to_perm_include, _raw_reference_path, _destination_folder):
        # First parameters
        self.n_refs = _n_refs
        self.refs_len = _refs_len
        self.bases_to_include = _bases_to_perm_include
        self.raw_reference_path = _raw_reference_path

        # We asume that VCFParser is excecuted first
        tmp_folder = os.path.join(_destination_folder, "Tmp")
        self.parsing_folder = os.path.join(tmp_folder, "Parsing")
        self.metadata_folder = os.path.join(tmp_folder, "Meta_data")

        assert os.path.isdir(tmp_folder)
        assert os.path.isdir(self.parsing_folder)
        assert os.path.isdir(self.metadata_folder)

        # Helper variables
        self.current_ref_has_lenght = False
        self.current_ref_data = {}

    def ProcessDescriptionLine(self, description_line):
        description_line = description_line[1:] # Remove ">" of the begining
        elements = description_line.split(" ")
        # Is mandatory that the first element will be the id
        self.current_ref_data["ID"] = elements[0]
        # Util description not always appear :(
        if len(elements) > 1:
            for element in elements[1:]:
                


    def StartReferenceProcessing(self):
        processed_ref_file_path = os.path.join(self.parsing_folder, "Reference.tmprlz")
        metadata_file_path = os.path.join(self.metadata_folder, "Reference.metarlz")
        
        is_first_time = True
        has_length = False
        description_line = ""
        with open(self.raw_reference_path, 'r') as raw_reference, open(processed_ref_file_path, 'wb') as processed_reference, open(metadata_file_path, 'wb') as metadata:
            # First we read the first line of the ref, which should be like >1 dna:chromosome chromosome:GRCh37:1:1:249250621:1\n
            description_line = raw_reference.readline()


            
        

