from VCFParser import VCFParser
import sys

if __name__ == "__main__":
    vcf_path = sys.argv[1]

    parser = VCFParser(vcf_path)

    parser.StartParsing()