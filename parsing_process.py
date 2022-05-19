from VCFParser import VCFParser
import sys

if __name__ == "__main__":
    vcf_path = sys.argv[1]
    # argv = [codigo, param1, param2, ...]

    if sys.argv[1] == "-n":
        # TODO: Ahora asumimos que todo correcto no mas
        parser = VCFParser(sys.argv[- int(sys.argv[2]):], True)
    else:
        # TODO
        exp = """
            Usage: parsing_process.py [OPTIONS]... [FILE]...\n
            For each VCF FILE, we parse given by the OPTIONS defined.\n
            Example: parsing_process.py -n 1 file.vcf\n\n
            Mandatory options in order:\n
            \t-n=NUM : Number of files, if NUM = 1, we assume all the genome is defined in the same file, else each file corresponds with one chromosome in order.\n
            \t-f=[STR]... :\n
            """
        print(exp)


    parser = VCFParser(vcf_path, True)

    parser.StartParsing()