import sys

if __name__ == "__main__":
    vcf_path = sys.argv[1]
    n_record = 0
    with open(vcf_path, 'r') as VCF:
        line = VCF.readline()

        while line[0:2] == "##":
            line = VCF.readline()
        
        header = line.split("\t")[9:]

        while line:
            n_record += 1
            line = VCF.readline()

        print("This VCF contains", n_record, "records.")