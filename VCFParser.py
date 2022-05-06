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

        self.path_fileParsed = VCF_path[:-4] + "_Parsed.txt"
        self.VCFParsed = None

        # Preference settings
        self.UnphasedAsPhased = True
        self.DiscardNotPASSRecords = True

        # Helper parameters
        self.p_metainfo_line = re.compile(r"^##")
        self.p_record_line = re.compile(r"^#")
        self.p_nucleotid_only = re.compile(r"[ACTGN]+")
        self.p_SVTYPE = re.compile(r"SVTYPE=[^;]+")
            # Case1 like AAA[<ctg>:2[, Case2 like AAA]<ctg>:2], Case3 like [<ctg>:2[AAA and Case4 like ]<ctg>:2]AAA
        self.p_bndCase1 = re.compile(r"[ACTGN]+\[<[0-9A-Za-z!#$%&+./:;?@^_|~-][0-9A-Za-z!#$%&*+./:;=?@^_|~-]*>:[0-9]+\[")
        self.p_bndCase2 = re.compile(r"[ACTGN]+\]<[0-9A-Za-z!#$%&+./:;?@^_|~-][0-9A-Za-z!#$%&*+./:;=?@^_|~-]*>:[0-9]+\]")
        self.p_bndCase3 = re.compile(r"\[<[0-9A-Za-z!#$%&+./:;?@^_|~-][0-9A-Za-z!#$%&*+./:;=?@^_|~-]*>:[0-9]+\[[ACTGN]+")
        self.p_bndCase4 = re.compile(r"\]<[0-9A-Za-z!#$%&+./:;?@^_|~-][0-9A-Za-z!#$%&*+./:;=?@^_|~-]*>:[0-9]+\][ACTGN]+")

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
        self.ID_samples = []
        self.n_samples = 0
        self.Lenght_Reference = 0

        # Variables for VCF record processing
        self.curr_Chrom = "X"
        self.curr_Pos = 0
        self.curr_ID = "X"
        self.curr_REF = "X"
        self.curr_Alt = "X"
        self.curr_Qual = "X"
        self.curr_Filter = "X"
        self.curr_Info = {}
        self.curr_Format = {}
        self.curr_AleleList = []
        self.curr_SVTYPE = "X"

    def ReferenceIndexTransform(self, index):
        return 2 * self.Lenght_Reference - index - 1

    def CalculateInvertedReference(self):
        # TODO: Recordar generar un reporte del META con el tema de la REF inv y meta_ReferenceValues
        len_reference = 0
        for ID in self.meta_ReferenceValues.keys():
            if ID:
                len_reference += int(self.meta_ReferenceValues.get(ID).get("length"))    
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
                dict_aux["ID"] = self.counter_contig # Set ID to a shorter internal value as new ID
                self.counter_contig += 1
                dict_aux.pop("ID")

                self.meta_ReferenceValues[ID] = dict_aux # = {'GL000224.1': {'ID': 1,'assembly': 'b37', 'length': '179693'}}
            
            if self.isDebugMode: print(self.meta_ReferenceValues)

            # TODO: The rest of the lines
            line = self.VCF.readline()[:-1]
        
        # This last line its supposed to be the header line
        self.ID_samples = line.split("\t")[9:]
        self.n_samples = len(self.ID_samples) 
        if self.isDebugMode: 
            print("Curr n_samples: ", self.n_samples)
            print("Curr IDs: ", self.ID_samples)
            print("Destination path: ", self.path_fileParsed)

        self.CalculateInvertedReference()


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
            return raw_AleleList.split("|") if self.UnphasedAsPhased else []
        else:
            return raw_AleleList.split("|")

    def GetSVTYPE(self, match_SVTYPE):
        start, end = match_SVTYPE.span()
        return match_SVTYPE.string[start + 7 : end]

    # def WritePhrase2(self):
    #     phrase = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t".format(self.phrase_INDV,
    #                                                         self.phrase_Chrom,
    #                                                         self.phrase_Alele,
    #                                                         self.phrase_Pos,
    #                                                         self.phrase_Len,
    #                                                         self.phrase_Edit,
    #                                                         self.phrase_PosEdit,
    #                                                         self.phrase_LenEdit)
    #     # TODO: Este sistema funciona?
    #     self.path_fileParsed.write(phrase)

    def WritePhrase(self, values_phrase):
        phrase = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t".format(values_phrase[0],
                                                            values_phrase[1],
                                                            values_phrase[2],
                                                            values_phrase[3],
                                                            values_phrase[4],
                                                            values_phrase[5],
                                                            values_phrase[6],
                                                            values_phrase[7])
        # TODO: Este sistema funciona?
        self.VCFParsed.write(phrase)

    def WriteInternalINV(self):
        pos_END = None
        pos_END = self.curr_Info.get("END")

        if pos_END:
            pos_END = int(pos_END) - 1
            # TODO: Si cualquier pos > lenRef, se interpreta la inversa
            self.phrase_Len = pos_END - self.curr_Pos
            self.phrase_Edit = self.meta_ReferenceValues.get(self.curr_Chrom).get("ID")
            self.phrase_PosEdit = self.ReferenceIndexTransform(pos_END)
            self.phrase_LenEdit = self.phrase_Len
        else:
            print("(!!) ERROR|INV: Se ha encontrado una inversion sin campo END, se descartará el edit.")
            return

        self.WritePhrase()

    def WriteDeletion(self):
        pos_END = int(self.curr_Info.get("END")) - 1 # Correction for 0 start

        self.phrase_Len = pos_END - self.phrase_Pos + 1 # I know +-1 is unnecesary, Its just for theorical coherence
        self.phrase_Edit = ""
        self.phrase_PosEdit = 0
        self.phrase_LenEdit = 0

        self.WritePhrase() # Done

    def AddToPhraseCache(self):
        tmp_phrase = [-1, # self.phrase_INDV to complete
                    self.phrase_Chrom,
                    -1, # self.phrase_Alele to complete
                    self.phrase_Pos,
                    self.phrase_Len,
                    self.phrase_Edit,
                    self.phrase_PosEdit,
                    self.phrase_LenEdit]
        self.phrase_Cache.append(tmp_phrase)

    def ProcessVariants(self):
        """
        Remains:
                self.phrase_INDV = "X"
                self.phrase_Alele = 0
                self.phrase_Edit = "X"
                self.phrase_PosEdit = 0
                self.phrase_LenEdit = 0
        """
        for alt in self.curr_AltList:
            if re.fullmatch(self.p_nucleotid_only, alt): # If its an explicit edit
                self.phrase_PosEdit = 0
                self.phrase_LenEdit = 0
                self.phrase_Edit = alt

                self.AddToPhraseCache() # Done
            else:
                print("Jaja perate")
                



    def ProcessRECORDS(self):
    
        raw_record = self.VCF.readline()

        while raw_record:
            record = raw_record[:-1].split('\t')

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
                self.phase_INDV = self.ID_samples[i]

                self.UpdateAlelesList(raw_AleleFullList[i])

                for j in range(len(self.curr_AleleList)): # Over each alele

                    if self.curr_AleleList[j] == 0: # If there's no change, we continue
                        continue

                    tmp_values_phrase = self.phrase_Cache[self.curr_AleleList[j] - 1]
                    tmp_values_phrase[0] = self.phrase_INDV
                    tmp_values_phrase[2] = j

                    self.WritePhrase(tmp_values_phrase)

                    # """
                    # Remains:
                    #     self.phrase_PosEdit = 0
                    #     self.phrase_LenEdit = 0
                    # """

                    # # Check the edit

                    # if re.fullmatch(self.p_nucleotid_only, self.phrase_Edit): # If its an explicit edit
                    #     self.phrase_PosEdit = 0
                    #     self.phrase_LenEdit = 0

                    #     self.WritePhrase() # Done
                    #     continue

                    # elif self.curr_Info.contains_key("SVTYPE"): # If its an external reference edit, we should check the SVTYPE
                    #     # TODO: Chequear si es necesario hacer una variable de clase
                    #     self.curr_SVTYPE = self.curr_Info.get("SVTYPE")
                    #     bool_HasEND = self.curr

                    #     # We need to check what kind if SVTYPE is
                    #         # No estoy segura de (self.curr_SVTYPE == "DEL" and self.phrase_Edit == "<DEL>")
                    #     if (self.phrase_Edit == "<DEL>" or self.phrase_Edit == "<CN0>"):
                    #         self.WriteDeletion()
                    #         continue
                        
                    #     # TODO: No estoy segura de si esto es cierto. VERIFICAR
                    #     # https://github.com/vcflib/vcflib/blob/master/src/Variant.cpp 359
                    #     elif (self.curr_SVTYPE == "INS"):
                    #         # TODO: Las inserciones e inversiones se ven aqui
                    #         # (!) Si X != REF => Interpretar dos edits
                    #         if re.match(self.p_bndCase1, self.phrase_Edit):
                    #             pass
                    #         elif re.match(self.p_bndCase2, self.phrase_Edit):
                    #             pass
                    #         elif re.match(self.p_bndCase3, self.phrase_Edit):
                    #             pass
                    #         elif re.match(self.p_bndCase4, self.phrase_Edit):
                    #             pass

                    #     elif (self.curr_SVTYPE == "INV"):
                    #         assert self.phrase_Edit == "<INV>"
                    #         self.WriteInternalINV()
                    #         continue

                    #     elif (self.curr_SVTYPE == "DUP"):
                    #         # TODO: Handle
                    #         pass

                    #     elif (self.curr_SVTYPE == "CNV"):
                    #         # TODO: Handle
                    #         pass

                    # else:
                    #     # TODO: Aqui caen todos los casos no manejados, crear un reporte de no soportados
                    #     pass
            raw_record = self.VCF.readline()

    def StartParsing(self):

        with open(self.path_fileParsed, 'w') as aux_VCFParsed:
            self.VCFParsed = aux_VCFParsed
            with open(self.path_file, 'r') as aux_VCF:
                
                self.VCF = aux_VCF
                # Collect VCF metainformation
                self.ProcessMETA()    
                # Interpretate edits
                self.ProcessRECORDS()        




