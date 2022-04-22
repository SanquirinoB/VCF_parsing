import re

class VCFParser():
    def __init__(self, VCF_path) -> None:
        # Reference to the VCF to be parsed
        self.file_path = VCF_path

        # Helper parameters
        self.p_metainfo_line = re.compile(r"^##")
        self.p_record_line = re.compile(r"^#")
        self.p_nucleotid_only = re.compile(r"[ACTGN]+")
        # TODO: Esto aun no funciona
        self.p_SVTYPE = re.compile(r"SVTYPE=.+[^;];?")

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
        self.curr_Info = "X"
        self.curr_Format = "X"
        self.curr_AleleList = []

    def WritePhrase(self):
        # TODO: Implementar escritura en archivo
        pass

    def StartParsing(self):
        with open(self.file_path, 'r') as VCF:

            # Collect VCF metainformation
            meta = VCF.readline()
            # Muy overkill while re.match(self.p_metainfo_line, meta):
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
                if (record[6] != "PASS"):
                    print("(!) WARNING|FILTER: Se ha descartado un registro por no cumplir con FILTER=PASS. Edit nro {}".format(i+1))
                    continue

                # Update internal variables
                self.curr_Chrom = int(record[0])
                self.curr_Pos = int(record[1])
                self.curr_ID = record[2]
                self.curr_REF = record[3]
                self.curr_AltList = record[4].split(',')
                self.curr_Qual = record[5]
                self.curr_Filter = record[6]
                self.curr_Info = record[7]
                self.curr_Format = record[8]
                raw_AleleFullList = record[9:]

                # Over each sample
                for i in range(n_samples):
                    self.phase_INDV = ID_samples[i]

                    # TODO: Por ahora asumimos que siempre hay phased
                    raw_AleleList = raw_AleleFullList[i].replace('/', '|')
                    self.curr_AleleList = [int(x) for x in raw_AleleList.split('|')]

                    # Over each alele
                    for j in range(len(self.curr_AleleList)):
                        self.phrase_Alele = self.curr_AleleList[j]

                        # If there's no change, we continue
                        if self.phrase_Alele == 0:
                            continue

                        self.phrase_Edit = self.curr_AltList[self.phrase_Alele - 1]

                        # Check the edit
                            # If its an explicit edit
                        if re.fullmatch(self.p_nucleotid_only, self.phrase_Edit):
                            self.phrase_PosEdit = 0
                            self.phrase_LenEdit = 0
                            # TODO: Donde me aseguro que todo tiene el valor que le corresponde?
                            self.WritePhrase()

                        else:
                            # (!) INFO: SVTYPE=BND
                            # <id_ALT>
                                # <DEL>
                                # <INS>
                                # <DUP>
                                # <INV>
                                # <CNV>

                            # (!) Si X != REF => Interpretar dos edits
                            # X]<>:p]
                            # X[<>:p[
                            # ]<>:p]X
                            # [<>:p[X

                raw_record = VCF.readline()


