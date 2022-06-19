#!/usr/bin/python
# -*- coding: utf-8 -*-

from VCFParser import VCFParser
import sys

def main():
    if sys.argv[3] == "-n":
        # TODO: Ahora asumimos que todo correcto no mas
        # (self, Destination_folder, VCF_path_list, MISS_AleleAlt_Value = 0, LeaveUnphasedAsPhased = True, 
        # DiscardNotPASSRecords = True, debug = False):
        parser = VCFParser(sys.argv[1], sys.argv[2], sys.argv[- int(sys.argv[4]):], 0, True, True, False)
        parser.StartParsing()
        return 0
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
        return 1

    

if __name__ == "__main__":
    sys.exit(main())

    