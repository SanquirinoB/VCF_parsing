from VCFParser import VCFParser
import sys

if __name__ == "__main__":
    # argv = [codigo, destino, -n, paths, ... , param2, ...]

    if sys.argv[2] == "-n":
        # TODO: Ahora asumimos que todo correcto no mas
        # (self, Destination_folder, VCF_path_list, MISS_AleleAlt_Value = 0, LeaveUnphasedAsPhased = True, 
        # DiscardNotPASSRecords = True, debug = False):
        parser = VCFParser(sys.argv[1], sys.argv[- int(sys.argv[2]):], 0, True, True, True)
    else:
        # TODO usar https://docs.python.org/3/library/argparse.html
        exp = """
            Usage: parsing_process.py [OPTIONS]... [FILE]...\n
            For each VCF FILE, we parse given by the OPTIONS defined.\n
            Example: parsing_process.py -n 1 file.vcf\n\n
            Mandatory options in order:\n
            \t-n=NUM : Number of files, if NUM = 1, we assume all the genome is defined in the same file, else each file corresponds with one chromosome in order.\n
            \t-f=[STR]... :\n
            """
        print(exp)
        exit(1)

    parser.StartParsing()