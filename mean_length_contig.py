import sys

if __name__ == "__main__":
    vcf_path = sys.argv[1]
    length, n = 0, 0
    with open(vcf_path, 'r') as VCF:
        line = VCF.readline()

        while line[:2] == "##":
            
            if line[:9] == "##contig=": # line = "##contig=<ID=GL000224.1,assembly=b37,length=179693>\n"

                for x in line[10:-1].split(","): # line[10:-1] = "ID=GL000224.1,assembly=b37,length=179693"  
                    pair = x.split("=")

                    if pair[0] == "length": # Necessary for invertion calculus
                        length += int(pair[1])
                        n += 1
            line = VCF.readline()

    print("Max pos over reference is arround", length / n, ".")