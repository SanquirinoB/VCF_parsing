import sys

if __name__ == "__main__":
    vcf_path = sys.argv[1]

    with open(vcf_path, 'r') as VCF:
        line = VCF.readline()

        while line[0:2] == "##":
            line = VCF.readline()
        
        header = line.split("\t")[9:]

        print("This VCF contains", len(header), "samples.")