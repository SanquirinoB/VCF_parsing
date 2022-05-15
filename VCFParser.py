#!/usr/bin/python
# -*- coding: utf-8 -*-

import re

class VCFParser():
    def __init__(self, VCF_path, debug):
        # Debug
        self.isDebugMode = debug
        # Reference to the VCF to be parsed
        self.path_file = VCF_path
        self.VCF = None
        self.MISSING = "."
        # TODO: Hacer configurable a gusto
        self.MISS_AleleAlt = 0

        self.path_fileParsed = VCF_path[:-4] + "_Parsed.txt"
        self.VCFParsed = None

        # Preference settings
        self.UnphasedAsPhased = True
        self.DiscardNotPASSRecords = True

        # Helper parameters
        self.p_metainfo_line = re.compile(r"^##")
        self.p_record_line = re.compile(r"^#")
        self.p_nucleotid_only = re.compile(r"[ACTGN]+")
        self.p_cnv_record = re.compile(r"<CN[1-9][0-9]*>") # => <CNi> where i >= 1
        self.p_SVTYPE = re.compile(r"SVTYPE=[^;]+")
            # Case1 like AAA[<ctg>:2[, Case2 like AAA]<ctg>:2], Case3 like [<ctg>:2[AAA and Case4 like ]<ctg>:2]AAA
        # self.p_bndCase1 = re.compile(r"[ACTGN]+\[<[0-9A-Za-z!#$%&+./:;?@^_|~-][0-9A-Za-z!#$%&*+./:;=?@^_|~-]*>:[0-9]+\[")
        # self.p_bndCase2 = re.compile(r"[ACTGN]+\]<[0-9A-Za-z!#$%&+./:;?@^_|~-][0-9A-Za-z!#$%&*+./:;=?@^_|~-]*>:[0-9]+\]")
        # self.p_bndCase3 = re.compile(r"\[<[0-9A-Za-z!#$%&+./:;?@^_|~-][0-9A-Za-z!#$%&*+./:;=?@^_|~-]*>:[0-9]+\[[ACTGN]+")
        # self.p_bndCase4 = re.compile(r"\]<[0-9A-Za-z!#$%&+./:;?@^_|~-][0-9A-Za-z!#$%&*+./:;=?@^_|~-]*>:[0-9]+\][ACTGN]+")

        # Phrase base structure
            # Everything related to indexes will be fixed for start at 0
            # for the first element
        self.phrase_INDV = "X"
        self.phrase_Chrom = 0
        self.phrase_Alele = 0
        self.phrase_Pos = 0
        self.phrase_Len = 0
        self.phrase_Edit = "X"
        self.phrase_PosEdit = 0
        self.phrase_LenEdit = 0

        self.phrase_Cache = []

        # Variables for VCF metadata processing
        self.meta_ReferenceValues = {}
        self.counter_contig = 0
        self.ID_samples = {}
        self.n_samples = 0
        self.Lenght_Reference = 0

        # Variables for VCF record processing
        self.curr_Chrom = "X"
        self.curr_Pos = 0
        self.curr_ID = "X"
        self.curr_REF = "X"
        self.curr_AltIndex = 0
        self.curr_Qual = "X"
        self.curr_Filter = "X"
        self.curr_Info = {}
        self.curr_Format = {}
        self.curr_AleleList = []
        self.curr_SVTYPE = "X"

    def ReferenceIndexTransform(self, index):
        return (2 * self.Lenght_Reference) - index - 1

    def CalculateInvertedReference(self):
        # TODO: Recordar generar un reporte del META con el tema de la REF inv y meta_ReferenceValues
        for ID in self.meta_ReferenceValues.keys():
            if ID:
                self.Lenght_Reference += int(self.meta_ReferenceValues.get(ID).get("length"))    
            else:
                exit(1)   

    def ProcessMETA(self):
        # The first line readed is the VCF Version
        # TODO: Querremos procesar esto? Quiza crear un assert de version
        line = self.VCF.readline()[:-1]
        pair, dict_aux = [], {}

        # The last line allowed will be just before header line
        while line[:2] == "##":

            if line[:9] == "##contig=": # line = "##contig=<ID=GL000224.1,assembly=b37,length=179693>\n"

                if not ("length" in line): # This value is not mandatory, so just in case
                    # TODO: Debemos recuperar el archivo de referencia y recuperar el largo
                    # de los strings
                    print("Oh no")

                for x in line[10:-1].split(","): # line[10:-1] = "ID=GL000224.1,assembly=b37,length=179693"  
                    pair = x.split("=")
                    dict_aux[pair[0]] = pair[1]

                ID = dict_aux.get("ID") # = {'ID': 'GL000224.1', 'assembly': 'b37', 'length': '179693'}
                dict_aux["internalID"] = self.counter_contig # Set ID to a shorter internal value as new ID
                self.counter_contig += 1

                self.meta_ReferenceValues[ID] = dict_aux # = {'GL000224.1': {'ID': 1,'assembly': 'b37', 'length': '179693'}}
            

            # TODO: The rest of the lines
            line = self.VCF.readline()[:-1]
        
        # This last line its supposed to be the header line
        tmp_ID_samples = line.split("\t")[9:]
        self.n_samples = len(tmp_ID_samples) 
        self.ID_samples = {}
        for i in range(self.n_samples):
            self.ID_samples[i] = tmp_ID_samples[i]
        
        self.CalculateInvertedReference()

        if self.isDebugMode: 
            print("---- RESUME META ----")
            print("\tCurr n_samples:", self.n_samples)
            print("\tCurr IDs:", self.ID_samples)
            print("\tDestination path:", self.path_fileParsed)
            print("\tMeta reference values:", self.meta_ReferenceValues)
            print("\tReference length:", self.Lenght_Reference)
        


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
        if raw_INFO == self.MISSING: return
        list_INFO = raw_INFO.split(";")
        dict_INFO = {}
        for infoPair in list_INFO:
            key, value = infoPair.split("=")

            if key == "END" or key == "SVLEN":
                value = [int(x) for x in value.split(",")]

            dict_INFO[key] = value

        self.curr_Info = dict_INFO

    def UpdateInternalValues(self, record):
        """
        Updates all internal variables which start with curr_*
        This information is needed for global record processing
        """
        self.curr_Chrom = record[0]
        self.curr_Pos = int(record[1])
        self.curr_ID = record[2]
        self.curr_REF = record[3]
        self.curr_AltList = record[4].split(',')
        if self.isDebugMode: print("Current AltList: ", self.curr_AltList)

        self.curr_Qual = record[5]
        self.curr_Filter = record[6]
        self.ProcessINFO(record[7])
        self.ProcessFORMAT(record[8])

    def UpdateGenericPhraseValues(self):
        self.phrase_Chrom = self.curr_Chrom
        self.phrase_Pos = self.curr_Pos
        self.phrase_Len = len(self.curr_REF)

    def UpdateAlelesList(self, INDVRecord):
        tmp_INDVRecord = INDVRecord.split(":")
        raw_AleleList = tmp_INDVRecord[self.curr_Format.get("GT")]

        if ("/" in raw_AleleList):
            raw_AleleList = raw_AleleList.replace("/", "|")
            self.curr_AleleList = raw_AleleList.split("|") if self.UnphasedAsPhased else []
        else:
            self.curr_AleleList = raw_AleleList.split("|")

        self.curr_AleleList = [int(x) if x != "." else self.MISS_AleleAlt for x in self.curr_AleleList]

        if self.isDebugMode: print("CurrAleleList is:", self.curr_AleleList)

    def GetSVTYPE(self, match_SVTYPE):
        start, end = match_SVTYPE.span()
        return match_SVTYPE.string[start + 7 : end]

    def WritePhrase(self, list_values_phrase):
        if self.isDebugMode: print("Phrases to be writed: ", len(list_values_phrase))

        for values_phrase in list_values_phrase:
            phrase = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(values_phrase[0],
                                                                values_phrase[1],
                                                                values_phrase[2],
                                                                values_phrase[3],
                                                                values_phrase[4],
                                                                values_phrase[5],
                                                                values_phrase[6],
                                                                values_phrase[7])
            
            if self.isDebugMode: print("Phrase to be writed: ", phrase)

            phrase = phrase.encode("utf-8")
            self.VCFParsed.write(phrase)
    
    def FindValidVariantLength(self):

        if "END" in self.curr_Info.keys():
            # END is 1-based, tecnically this return should be (END - 1) - POS + 1
            return self.curr_Info.get("END")[self.curr_AltIndex - 1] - self.curr_Pos
        elif "SVLEN" in self.curr_Info.keys():
            return self.curr_Info.get("SVLEN")[self.curr_AltIndex - 1] 
        else:
            # TODO: Handle Error
            print("Oh no, error FindValidVariantLength")
            exit(1)
        
    def GenerateDeletionPhraseCache(self):
        edit_length = self.FindValidVariantLength() # If it doesnt work, the exception will be thrown in this function
        self.phrase_Len = edit_length
        self.phrase_Edit = ""
        self.phrase_PosEdit = 0
        self.phrase_LenEdit = 0

        self.AddToPhraseCache() # Done

    def GenerateInversionPhraseCache(self):
        edit_length = self.FindValidVariantLength() # If it doesnt work, the exception will be thrown in this function
        self.phrase_Len = edit_length
        self.phrase_Edit = self.meta_ReferenceValues.get(self.curr_Chrom).get("ID")
        self.phrase_PosEdit = self.ReferenceIndexTransform(self.curr_Pos + edit_length)
        self.phrase_LenEdit = self.phrase_Len

        self.AddToPhraseCache() # Done

    def GenerateDuplicationPhraseCache(self):
        edit_length = self.FindValidVariantLength() # If it doesnt work, the exception will be thrown in this function
        n_copy = int(self.curr_Alt[3:-1])
        
        self.phrase_Len = 0
        self.phrase_Edit = self.meta_ReferenceValues.get(self.curr_Chrom).get("ID")
        self.phrase_PosEdit = self.curr_Pos
        self.phrase_LenEdit = edit_length

        list_tmp_phrase = [self.CreateCustomPhrase() for _ in range(n_copy)]
        self.CustomAddToPhraseCache(list_tmp_phrase)

    def AddToPhraseCache(self):
        tmp_phrase = [-1, # self.phrase_INDV to complete
                    self.phrase_Chrom,
                    -1, # self.phrase_Alele to complete
                    self.phrase_Pos,
                    self.phrase_Len,
                    self.phrase_Edit,
                    self.phrase_PosEdit,
                    self.phrase_LenEdit]
        self.CustomAddToPhraseCache([tmp_phrase])

    def CustomAddToPhraseCache(self, tmp_phrase):
        self.phrase_Cache.append(tmp_phrase)
    
    def CreateCustomPhrase(self, Chrom = None, Pos = None, Len = None,
                            Edit = None, PosEdit = None, LenEdit = None):
        tmp_phrase = [-1, # self.phrase_INDV to complete
                    Chrom if Chrom else self.phrase_Chrom,
                    -1, # self.phrase_Alele to complete
                    Pos if Pos else self.phrase_Pos,
                    Len if Len else self.phrase_Len,
                    Edit if Edit else self.phrase_Edit,
                    PosEdit if PosEdit else self.phrase_PosEdit,
                    LenEdit if LenEdit else self.phrase_LenEdit]
        return tmp_phrase

    def ProcessVariants(self):

        for alt in self.curr_AltList:
            self.curr_Alt = alt
            if re.fullmatch(self.p_nucleotid_only, alt): # If its an explicit edit
                self.phrase_PosEdit = 0
                self.phrase_LenEdit = 0
                self.phrase_Edit = alt

                self.AddToPhraseCache() # Done
            elif "SVTYPE" in self.curr_Info.keys(): # If its an external reference edit, we should check the SVTYPE
                self.curr_SVTYPE = self.curr_Info.get("SVTYPE")

                if self.curr_SVTYPE == "DEL" or alt == "<CN0>":
                    self.GenerateDeletionPhraseCache()
                    
                elif self.curr_SVTYPE == "INV" and alt == "<INV>":
                    self.GenerateInversionPhraseCache()

                elif self.curr_SVTYPE == "DUP" and re.fullmatch(self.p_cnv_record, alt):
                    self.GenerateDuplicationPhraseCache()

            else:
                print("Jaja perate")
                exit(1)
    
    def CleanUpData(self):
        self.phrase_Cache = []

    def ProcessRECORDS(self):
    
        raw_record = self.VCF.readline()

        while raw_record:
            record = raw_record[:-1].split('\t')

            self.CleanUpData()

            if self.isDebugMode: print("Current record: ", record)


            # Filter check
            if (self.DiscardNotPASSRecords and record[5] != "PASS"):
                print("(!) WARNING|FILTER: Se ha descartado un registro por no cumplir con FILTER=PASS. Edit nro {}".format(-1))
                raw_record = self.VCF.readline()
                continue

            self.UpdateInternalValues(record)

            raw_AleleFullList = record[9:]

            self.UpdateGenericPhraseValues()

            self.ProcessVariants()

            # Over each sample
            for i in range(self.n_samples):
                if self.isDebugMode: print("For sample ", self.ID_samples.get(i))
                # Set internal id
                self.phrase_INDV = i

                self.UpdateAlelesList(raw_AleleFullList[i])

                for j in range(len(self.curr_AleleList)): # Over each alele
                    if self.isDebugMode: print("For alele ", j, "with variant index", self.curr_AleleList[j] - 1)
                    
                    if self.curr_AleleList[j] == 0: # If there's no change, we continue
                        continue
                    

                    self.curr_AltIndex = self.curr_AleleList[j] - 1
                    tmp_values_phrase = self.phrase_Cache[self.curr_AltIndex]

                    if self.isDebugMode: print("Phrases to be writed:", tmp_values_phrase)

                    for values_phrase in tmp_values_phrase:
                        values_phrase[0] = self.phrase_INDV
                        values_phrase[2] = j

                    self.WritePhrase(tmp_values_phrase)

            raw_record = self.VCF.readline()

    def StartParsing(self):

        with open(self.path_file, mode="r") as aux_VCF, open(self.path_fileParsed, mode="wb") as aux_VCFParsed:
            
            self.VCF = aux_VCF
            self.VCFParsed = aux_VCFParsed
            # Collect VCF metainformation
            self.ProcessMETA()    
            # Interpretate edits
            self.ProcessRECORDS()        
