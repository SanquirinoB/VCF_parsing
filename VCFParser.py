import re
from unittest import case

class VCFParser():
    def __init__(self, VCF_path) -> None:
        # Reference to the VCF to be parsed
        self.file_path = VCF_path

        # Preference settings
        self.UnphasedAsPhased = True
        self.DiscardNotPASSRecords = True

        # Helper parameters
        self.p_metainfo_line = re.compile(r"^##")
        self.p_record_line = re.compile(r"^#")
        self.p_nucleotid_only = re.compile(r"[ACTGN]+")
        self.p_SVTYPE = re.compile(r"SVTYPE=[^;]+")

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

        # Variables for VCF metadata processing

        # Variables for VCF record processing
        self.curr_Chrom = 0
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
        self.curr_Chrom = int(record[0])
        self.curr_Pos = int(record[1])
        self.curr_ID = record[2]
        self.curr_REF = record[3]
        self.curr_AltList = record[4].split(',')
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

    def WritePhrase(self):
        # TODO: Implementar escritura en archivo
        pass

    def StartParsing(self):
        with open(self.file_path, 'r') as VCF:

            # Collect VCF metainformation

            meta = VCF.readline()

            while meta[0:2] == "##":
                # TODO: Save data
                meta = VCF.readline()

            header = meta

            ID_samples = header[9:]
            n_samples = len(ID_samples)

            raw_record = VCF.readline()

            while raw_record:
                record = raw_record.split('\t')

                # Filter check
                if (self.DiscardNotPASSRecords and record[6] != "PASS"):
                    print("(!) WARNING|FILTER: Se ha descartado un registro por no cumplir con FILTER=PASS. Edit nro {}".format(i+1))
                    continue

                self.UpdateInternalValues(record)
                raw_AleleFullList = record[9:]

                self.UpdateGenericPhraseValues()
                """
                Values setted at this point:
                    self.phrase_Chrom = self.curr_Chrom
                    self.phrase_Pos = self.curr_Pos
                    self.phrase_Len = len(self.curr_REF)

                Remains:
                    self.phrase_INDV = "X"
                    self.phrase_Alele = 0
                    self.phrase_Edit = "X"
                    self.phrase_PosEdit = 0
                    self.phrase_LenEdit = 0
                """

                # Over each sample
                for i in range(n_samples):
                    self.phase_INDV = ID_samples[i]
                    """
                    Remains:
                        self.phrase_Alele = 0
                        self.phrase_Edit = "X"
                        self.phrase_PosEdit = 0
                        self.phrase_LenEdit = 0
                    """

                    self.UpdateAlelesList(raw_AleleFullList[i])

                    for j in range(len(self.curr_AleleList)): # Over each alele

                        if self.curr_AleleList[j] == 0: # If there's no change, we continue
                            continue
                        
                        # Save Values
                        self.phrase_Alele = j
                        self.phrase_Edit = self.curr_AltList[self.phrase_Alele]
                        """
                        Remains:
                            self.phrase_PosEdit = 0
                            self.phrase_LenEdit = 0
                        """

                        # Check the edit

                        if re.fullmatch(self.p_nucleotid_only, self.phrase_Edit): # If its an explicit edit
                            self.phrase_PosEdit = 0
                            self.phrase_LenEdit = 0

                            self.WritePhrase() # Done
                            continue

                        elif self.curr_Info.contains_key("SVTYPE"): # If its an external reference edit, we should check the SVTYPE
                            # TODO: Chequear si es necesario hacer una variable de clase
                            self.curr_SVTYPE = self.curr_Info.get("SVTYPE")

                            # We need to check what kind if SVTYPE is
                            if (self.curr_SVTYPE == "DEL" and self.phrase_Edit == "<DEL>"):
                                pos_END = int(self.curr_Info.get("END")) - 1 # Correction for 0 start

                                self.phrase_Len = pos_END - self.phrase_Pos + 1 # I know +-1 is unnecesary, Its just for theorical coherence
                                self.phrase_Edit = ""
                                self.phrase_PosEdit = 0
                                self.phrase_LenEdit = 0

                                self.WritePhrase() # Done
                                continue

                            elif (self.curr_SVTYPE == "INS"):
                                # TODO: Handle
                                # (!) Si X != REF => Interpretar dos edits
                                    # X]<>:p]
                                    # X[<>:p[
                                    # ]<>:p]X
                                    # [<>:p[X
                                    # <wea>
                                pass
                            elif (self.curr_SVTYPE == "DUP"):
                                # TODO: Handle
                                pass
                            elif (self.curr_SVTYPE == "INV"):
                                # TODO: Handle
                                pass
                            elif (self.curr_SVTYPE == "CNV"):
                                # TODO: Handle
                                pass
                            elif (self.curr_SVTYPE == "BND"):
                                # TODO: Handle
                                pass

                        else:
                            # TODO: Aqui caen todos los casos no manejados, crear un reporte de no soportados
                            pass

                raw_record = VCF.readline()


