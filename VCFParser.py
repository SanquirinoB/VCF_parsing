#!/usr/bin/python
# -*- coding: utf-8 -*-

import re
import os
import shutil
from collections import OrderedDict
from Structures import Phrase, MetaInfo
from ReferenceProcessor import ReferenceProcessor


class VCFParser():
    def __init__(self, Destination_folder, Reference_path, VCF_path_list, N_chromosomes="24", MISS_AleleAlt_Value=0, LeaveUnphasedAsPhased=True,
                 DiscardNotPASSRecords=True, debug=False):
        # Debug
        self.isDebugMode = debug
        # Reference to the VCF to be parsed
        self.path_destination = Destination_folder
        self.path_file_list = VCF_path_list
        self.VCF = None
        self.MISSING = "."

        self.MISS_AleleAlt = MISS_AleleAlt_Value

        self.n_droppedRecords = 0
        self.size_VCFFiles_processed = 0
        self.size_ParsingFile_generated = 0
        self.VCFParsed = None

        # Preference settings
        self.UnphasedAsPhased = LeaveUnphasedAsPhased
        self.DiscardNotPASSRecords = DiscardNotPASSRecords

        # Helper parameters
        self.p_nucleotid_only = re.compile(r"[ACTGN]+")
        self.p_cnv_record = re.compile(
            r"<CN[1-9][0-9]*>")  # => <CNi> where i >= 1
        self.is_valid_record = True

        # Phrase base structure
        # Everything related to indexes will be fixed for start at 0
        # for the first element
        self.phrase_INDV = 0
        self.phrase_Chrom = 0
        # self.phrase_Alele = 0
        self.phrase_Pos = 0
        self.phrase_Len = 0
        self.phrase_Edit = 0
        self.phrase_PosEdit = 0
        self.phrase_LenEdit = 0

        self.n_print = 0

        self.phrase_Cache = []
        self.phrase_struct = Phrase()
        self.meta_struct = MetaInfo()

        # Variables for VCF metadata processing
        self.meta_ReferenceValues = OrderedDict()
        self.n_chromosomes = N_chromosomes
        self.ID_samples = {}
        self.n_samples = 0
        self.Length_Reference = 0
        self.n_phrases = 0
        self.reference_processor = ReferenceProcessor(
            Reference_path, Destination_folder)

        # Variables for VCF record processing
        self.curr_Chrom = ""
        self.curr_Pos = 0
        self.curr_REF = "X"
        self.curr_AltIndex = 0
        self.curr_Info = {}
        self.curr_Format = {}
        self.curr_AleleList = []
        self.curr_SVTYPE = "X"

        self.alphabet_replace = [["A", "1"], [
            "C", "2"], ["T", "3"], ["G", "4"], ["N", "5"]]

    def ReferenceIndexTransform(self, index):
        return (2 * self.Length_Reference) - index - 1

    def ProcessMETA(self, keep_meta=True):
        # The first line readed is the VCF Version
        # TODO: Querremos procesar esto? Quiza crear un assert de version
        line = self.VCF.readline()[:-1]

        # The last line allowed will be just before header line
        while line[:2] == "##":
            # All the info we actually need is provided by the Reference, but keep this process just in case of being necessary
            if keep_meta:
                pass
            line = self.VCF.readline()

        self.Length_Reference, self.meta_ReferenceValues = self.reference_processor.GetReferenceData()

        # This last line its supposed to be the header line
        tmp_ID_samples = line.split("\t")[9:]
        self.n_samples = len(tmp_ID_samples)
        # print(self.n_samples)
        self.ID_samples = {}
        for i in range(self.n_samples):
            self.ID_samples[i] = tmp_ID_samples[i]

    def ProcessFORMAT(self, raw_FORMAT):
        """
        Creates the dictionary of indexes for the sample data interpretation
        :param raw_FORMAT: String usually with the XX:XX:XX:xx format, where XX cuold be anything.
        """
        list_Format = raw_FORMAT.split(":")
        dict_FormatIndex = {}
        # So when we want to recover GT, just need to ask FORMATRecord.split(:)[?]
        # where ? is dict_FormatIndex.get("GT")
        for i in range(len(list_Format)):
            dict_FormatIndex[list_Format[i]] = i

        self.curr_Format = dict_FormatIndex

    def ProcessINFO(self, raw_INFO):
        if raw_INFO == self.MISSING:
            return
        list_INFO = raw_INFO.split(";")
        dict_INFO = {}
        key, value = "X", 0
        for infoPair in list_INFO:

            if "=" in infoPair:
                key, value = infoPair.split("=")
            else:  # Is flag
                key, value = infoPair, "True"

            if key == "END":
                value = [int(x) - 1 for x in value.split(",")]
            elif key == "SVLEN":
                value = [int(x) for x in value.split(",")]

            dict_INFO[key] = value

        self.curr_Info = dict_INFO

    def UpdateInternalValues(self, record):
        """
        Updates all internal variables which start with curr_*
        This information is needed for global record processing
        """
        self.curr_Chrom = record[0]
        self.curr_Pos = int(record[1]) - 1
        self.curr_REF = record[3]
        self.curr_AltList = record[4].split(',')

        # if self.isDebugMode:
        #     print("Current AltList: ", self.curr_AltList)

        self.ProcessINFO(record[7])
        self.ProcessFORMAT(record[8])

    def UpdateGenericPhraseValues(self):
        curr_Ref_data = self.meta_ReferenceValues.get(self.curr_Chrom)
        self.phrase_Chrom = int(self.curr_Chrom) - 1  # 0-index
        self.phrase_Pos = self.curr_Pos
        self.phrase_Len = len(self.curr_REF)

    def UpdateAlelesList(self, INDVRecord):
        tmp_INDVRecord = INDVRecord.split(":")
        raw_AleleList = tmp_INDVRecord[self.curr_Format.get("GT")]

        if ("/" in raw_AleleList):
            raw_AleleList = raw_AleleList.replace("/", "|")
            self.curr_AleleList = raw_AleleList.split(
                "|") if self.UnphasedAsPhased else []
        else:
            self.curr_AleleList = raw_AleleList.split("|")

        self.curr_AleleList = [
            int(x) if x.isnumeric() else self.MISS_AleleAlt for x in self.curr_AleleList]
        # if self.isDebugMode:
        #     print("CurrAleleList is:", self.curr_AleleList)

    def WritePhrase(self, list_values_phrase):
        # if True: print("Phrases to be writed: ", len(list_values_phrase))
        for values_phrase in list_values_phrase:
            self.phrase_struct.m_indv = values_phrase[0]
            self.phrase_struct.m_chrom = values_phrase[1]
            self.phrase_struct.m_alele = values_phrase[2]
            self.phrase_struct.m_pos = values_phrase[3]
            self.phrase_struct.m_len = values_phrase[4]
            self.phrase_struct.m_pos_e = values_phrase[6]
            self.phrase_struct.m_len_e = values_phrase[7]

            # if True: print("Phrase to be writed: ", values_phrase)

            #phrase = phrase.encode("ascii")
            self.VCFParsed.write(bytearray(self.phrase_struct))
            self.n_phrases += 1

            #self.n_print -= 1

    def FindValidVariantLength(self):

        if "END" in self.curr_Info.keys():
            # END is 1-based, tecnically this return should be (END - 1) - POS + 1
            return self.curr_Info.get("END")[self.curr_AltIndex - 1] - self.curr_Pos
        elif "SVLEN" in self.curr_Info.keys():
            return self.curr_Info.get("SVLEN")[self.curr_AltIndex - 1]
        else:
            print(
                "[FindValidVariantLength] ERROR: Valid end value (END/SVLEN) not found. Variant can't be processed.")
            exit(1)

    def CreateExplicitPhrase(self, edit):
        curr_ref_data = self.meta_ReferenceValues.get(edit)
        self.phrase_Edit = curr_ref_data.get("internal_ID")
        self.phrase_PosEdit = curr_ref_data.get("relPosRef")
        self.phrase_LenEdit = len(edit)

        return self.CreateCustomPhrase()  # Done

    def GenerateExplicitEditPhraseCacheFull(self, complete_edit):
        tmp_phrases_list = []

        while (len(complete_edit) > 4):
            tmp_phrases_list.append(
                self.CreateExplicitPhrase(complete_edit[:4]))
            complete_edit = complete_edit[4:]

        if (len(complete_edit) > 0):
            tmp_phrases_list.append(self.CreateExplicitPhrase(complete_edit))

        self.CustomAddToPhraseCache(tmp_phrases_list)

    def GenerateDeletionPhraseCache(self):
        # If it doesnt work, the exception will be thrown in this function
        edit_length = self.FindValidVariantLength()
        self.phrase_Len = edit_length
        self.phrase_Edit = 0
        self.phrase_PosEdit = 0
        self.phrase_LenEdit = 0

        self.AddToPhraseCache()  # Done

    def GenerateInversionPhraseCache(self):
        curr_Ref_data = self.meta_ReferenceValues.get(self.curr_Chrom)

        # If it doesnt work, the exception will be thrown in this function
        edit_length = self.FindValidVariantLength()
        self.phrase_Len = edit_length
        self.phrase_Edit = curr_Ref_data.get("internal_ID")
        # TODO no deberiamos trabajar con pos relativa, se arregla dsps
        self.phrase_PosEdit = self.ReferenceIndexTransform(
            curr_Ref_data.get("relPosRef") + self.curr_Pos + edit_length)
        self.phrase_LenEdit = self.phrase_Len

        self.AddToPhraseCache()  # Done

    def GenerateDuplicationPhraseCache(self):
        # If it doesnt work, the exception will be thrown in this function
        edit_length = self.FindValidVariantLength()
        n_copy = int(self.curr_Alt[3:-1])

        self.phrase_Edit = self.meta_ReferenceValues.get(
            self.curr_Chrom).get("internal_ID")
        self.phrase_PosEdit = self.curr_Pos
        self.phrase_LenEdit = edit_length

        tmp_phrases_list = [self.CreateCustomPhrase() for _ in range(n_copy)]
        self.CustomAddToPhraseCache(tmp_phrases_list)

    def AddToPhraseCache(self):
        tmp_phrase = [-1,  # self.phrase_INDV to complete
                      self.phrase_Chrom,
                      -1,  # self.phrase_Alele to complete
                      self.phrase_Pos,
                      self.phrase_Len,
                      self.phrase_Edit,
                      self.phrase_PosEdit,
                      self.phrase_LenEdit]
        self.CustomAddToPhraseCache([tmp_phrase])

    def CustomAddToPhraseCache(self, tmp_phrases_list):
        self.phrase_Cache.append(tmp_phrases_list)

    def CreateCustomPhrase(self, Chrom=None, Pos=None, Len=None,
                           Edit=None, PosEdit=None, LenEdit=None):
        tmp_phrase = [-1,  # self.phrase_INDV to complete
                      Chrom if Chrom else self.phrase_Chrom,
                      -1,  # self.phrase_Alele to complete
                      Pos if Pos else self.phrase_Pos,
                      Len if Len else self.phrase_Len,
                      Edit if Edit else self.phrase_Edit,
                      PosEdit if PosEdit else self.phrase_PosEdit,
                      LenEdit if LenEdit else self.phrase_LenEdit]

        return tmp_phrase

    def HasValidSTYPE(self):
        if ("SVTYPE" in self.curr_Info.keys()):
            self.curr_SVTYPE = self.curr_Info.get("SVTYPE")

            return (self.curr_SVTYPE == "DEL" or
                    self.curr_SVTYPE == "INV" or
                    self.curr_SVTYPE == "CNV" or
                    self.curr_SVTYPE == "DUP")
        else:
            return False

    def ProcessVariants(self):
        self.curr_AltIndex = 0
        for alt in self.curr_AltList:
            self.curr_Alt = alt

            if re.fullmatch(self.p_nucleotid_only, alt):  # If its an explicit edit
                self.GenerateExplicitEditPhraseCacheFull(alt)

            elif self.HasValidSTYPE():  # If its an external reference edit, we should check the SVTYPE

                if self.curr_SVTYPE == "DEL" or alt == "<CN0>":
                    self.GenerateDeletionPhraseCache()

                # TODO: Not supported for now
                elif self.curr_SVTYPE == "INV" and alt == "<INV>":
                    # self.GenerateInversionPhraseCache()
                    self.n_droppedRecords += 1
                    self.is_valid_record = False

                # SVTYPE tipically is CNV or DUP
                elif re.fullmatch(self.p_cnv_record, alt):
                    self.GenerateDuplicationPhraseCache()

                else:
                    self.n_droppedRecords += 1
                    self.is_valid_record = False

            else:
                if self.isDebugMode:
                    print("(!) Edit no canonico descartado.")
                self.n_droppedRecords += 1
                self.is_valid_record = False
            self.curr_AltIndex += 1

        if self.isDebugMode:
            print("Edits obtained: ", len(self.phrase_Cache),
                  "\nPhrases: ", self.phrase_Cache)

    def CleanUpData(self):
        self.phrase_Cache = []

    def ProcessRECORDS(self):

        raw_record = self.VCF.readline()

        while raw_record:
            record = raw_record[:-1].split('\t')

            # Fix temporal por samples
            if(len(record) < 10): break

            self.CleanUpData()

            # Filter check
            if (self.DiscardNotPASSRecords and record[6] != "PASS"):
                print(
                    "(!) WARNING|FILTER: Se ha descartado un registro por no cumplir con FILTER=PASS. Edit nro {}".format(-1))
                raw_record = self.VCF.readline()
                continue

            self.UpdateInternalValues(record)

            raw_AleleFullList = record[9:]

            self.UpdateGenericPhraseValues()

            self.ProcessVariants()

            if self.is_valid_record:
                # Over each sample
                for i in range(self.n_samples):
                    self.phrase_INDV = i

                    self.UpdateAlelesList(raw_AleleFullList[i])

                    for j in range(len(self.curr_AleleList)):  # Over each alele
                        #if self.isDebugMode: print("For alele ", j, "with variant index", self.curr_AleleList[j] - 1)

                        # If there's no change, we continue
                        if self.curr_AleleList[j] == 0:
                            continue

                        self.curr_AltIndex = self.curr_AleleList[j] - 1

                        tmp_values_phrase = self.phrase_Cache[self.curr_AltIndex]

                        #if self.isDebugMode: print("Phrases to be writed:", tmp_values_phrase)

                        for values_phrase in tmp_values_phrase:
                            values_phrase[0] = self.phrase_INDV
                            values_phrase[2] = j

                        self.WritePhrase(tmp_values_phrase)
            else:
                self.is_valid_record = True

            raw_record = self.VCF.readline()

    def GenerateMetaInfoResume(self):
        # Number of phrases (int)
        self.meta_struct.m_nPhrases = self.n_phrases
        self.TMPRLZ.write(bytearray(self.meta_struct))
        # Number of chromosomes (int)
        self.meta_struct.m_nPhrases = int(self.n_chromosomes)
        self.TMPRLZ.write(bytearray(self.meta_struct))

    def GenerateRLZResume(self):
        # Samples names and theirs IDs (IID, ID)
        # ID is related to the number line 1-indexed
        for key in self.ID_samples.keys():
            values = self.ID_samples[key]
            aux_line = "{}\r\n".format(values)
            self.TMPRLZ.write(aux_line.encode("utf8"))

    def ReportEndProcess(self):
        print("[RLZ] Resume parsing process:")
        print("[RLZ]\tNumber of samples identified:", self.n_samples)
        print("[RLZ]\tNumber of dropped records:", self.n_droppedRecords)
        print("[RLZ]\tNumber of valid edits captured:", self.n_phrases)
        print("[RLZ]\tFile size processed:",
              self.size_VCFFiles_processed / (1024 * 1024), "MB")
        print("[RLZ]\tSize of parsed file generated:",
              self.size_ParsingFile_generated / (1024 * 1024), "MB")

        if self.isDebugMode:
            print("\tCurr IDs:", self.ID_samples)
            print("\tMeta reference values:", self.meta_ReferenceValues)
            print("\tReference length:", self.Length_Reference)

    def StartParsing(self):
        # SUPUESTOS: Recibimos VCF separados por cromosoma

        tmp_folder = os.path.join(self.path_destination, "Tmp")
        parsing_folder = os.path.join(tmp_folder, "Parsing")
        metadata_folder = os.path.join(tmp_folder, "Meta_data")

        if os.path.isdir(tmp_folder):
            shutil.rmtree(tmp_folder)
        os.mkdir(tmp_folder)

        # TODO: Creo que estos dos ifs son innecesarios
        if os.path.isdir(parsing_folder):
            shutil.rmtree(parsing_folder)
        os.mkdir(parsing_folder)

        if os.path.isdir(metadata_folder):
            shutil.rmtree(metadata_folder)
        os.mkdir(metadata_folder)

        file_name = "phrases"
        path_fileParsed = os.path.join(
            parsing_folder, file_name + ".tmprlz")

        is_first_meta = True
        self.VCFParsed = open(path_fileParsed, mode="wb")

        self.reference_processor.StartReferenceProcessing()

        for path_file in self.path_file_list:
            file_stats = os.stat(path_file)
            self.size_VCFFiles_processed += file_stats.st_size
            # Get name and remove .vcf
            with open(path_file, mode="r") as aux_VCF:

                # Save pointers
                self.VCF = aux_VCF

                # Collect VCF metainformation
                self.ProcessMETA(keep_meta=is_first_meta)

                # First create a dummy start
                if (is_first_meta) : self.WritePhrase([[0, 0, 0, 0, 0, 0, 0, 0]])
                is_first_meta = False

                # Interpretate edits
                self.ProcessRECORDS()
            print("[RLZ]    ", path_file, "parsed!")

        # First create a dummy end
        self.WritePhrase([[self.n_samples - 1,
                           int(self.n_chromosomes) - 1,
                           1,
                           self.meta_ReferenceValues.get(
                               self.n_chromosomes).get("length") - 1,
                           0,
                           0,
                           0,
                           0]])

        self.VCFParsed.close()
        file_stats = os.stat(path_fileParsed)
        self.size_ParsingFile_generated = file_stats.st_size
        self.ReportEndProcess()

        path_ID_info = os.path.join(metadata_folder, "ID_info.metarlz")
        path_meta_info = os.path.join(metadata_folder, "Meta_info.metarlz")

        with open(path_meta_info, mode="wb") as aux_TMPRLZ:
            self.TMPRLZ = aux_TMPRLZ
            self.GenerateMetaInfoResume()

        with open(path_ID_info, mode="wb") as aux_TMPRLZ:
            self.TMPRLZ = aux_TMPRLZ
            self.GenerateRLZResume()
