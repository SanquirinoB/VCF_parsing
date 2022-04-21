import re

vcf_path = ''
p_nucleotid_only = re.compile(r"[ACTGN]+")

with open(vcf_path, 'r') as VCF:
    meta = VCF.readline()
    # Collect all meta information
    while meta[0:1] != "##":
        meta = VCF.readline()
    
    # Header line
    header = VCF.readline()
    id_samples = header[9:]
    n_samples = len(id_samples)
    # (I_i, C_i, A_i, p, l, E, p_i, l)
    phrase_aux = [0, 0, 0, 0, 0, 0, 0, 0]

    # Collect entries
    raw_data = VCF.readline()
    while raw_data:
        data = raw_data.split('\t')
        chrom = int(data[0])
        pos = int(data[1])
        id = data[2]
        ref = data[3]
        alt = data[4].split(',')
        qual = data[5]
        filter = data[6]
        info = data[7]
        format = data[8]
        id_haplo = data[9:]

        # (I_i, C_i, A_i, p, l, E, p_i, l)
        phrase_aux = ("X", chrom, -1, pos, len(ref), "X", 0, 0)

        # Over samples
        for i in range(1, n_samples + 1):
            # Set id of sample
            phrase_aux[0] = i
            id_haplo_aux = id_haplo[i].replace('/', '|')
            aleles = [int(x) for x in id_haplo_aux.split('|')]

            # Over aleles
            for alele_i in range(len(aleles)):
                alele_alt = aleles[alele_i]
                # If the alele uses the reference
                if alele_alt == 0:
                    continue

                # Save the alele edit applied
                phrase_aux[2] = alele_i

                # If the alele uses a variant
                    # Fix value to a valid array index
                alele_alt -=1
                variant = alt[alele_alt]

                    # Check if the variant is an explicit one
                if re.fullmatch(p_nucleotid_only, variant):
                    phrase_aux[5] = variant
                    continue
                    # If its not
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

        raw_data = VCF.readline()
