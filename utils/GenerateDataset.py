import sys
import os
from VCFUtils import VCFUtils

if __name__ == "__main__":
    destPath = sys.argv[1]
    util = VCFUtils()
    originalNames = []

    # with open(vcfList, 'r') as VCFList:
    #     line = VCFList.readline()
    #     print(line)
    #     while line:
    #         originalNames.append(line[:-1])
    #         line = VCFList.readline()

    # print(originalNames)
    # originalSizes = util.GetSizes(originalNames)
    originalSizes = [#["/d2/fernanda/vcf_base/c18.vcf", 2267185],
    #                  ["/d2/fernanda/vcf_base/c21.vcf", 1105538],
                     ["/d2/fernanda/vcf_base/c3.vcf", 5832277]
                    #  ["/d2/fernanda/vcf_base/c10.vcf", 3992219],
                    #  ["/d2/fernanda/vcf_base/c6.vcf", 5024119],
                    #  ["/d2/fernanda/vcf_base/c9.vcf", 3560687],
                    #  ["/d2/fernanda/vcf_base/c12.vcf", 3868428],
                    #  ["/d2/fernanda/vcf_base/c20.vcf", 1812841],
                    #  ["/d2/fernanda/vcf_base/c7.vcf", 4716715],
                    #  ["/d2/fernanda/vcf_base/c16.vcf", 2697949],
                    #  ["/d2/fernanda/vcf_base/c17.vcf", 2329288],
                    #  ["/d2/fernanda/vcf_base/c15.vcf", 2424689],
                    #  ["/d2/fernanda/vcf_base/c22.vcf", 1103547],
                    #  ["/d2/fernanda/vcf_base/c8.vcf", 4597105],
                    #  ["/d2/fernanda/vcf_base/c4.vcf", 5732585],
                    #  ["/d2/fernanda/vcf_base/c5.vcf", 5265763],
                    #  ["/d2/fernanda/vcf_base/c1.vcf", 6468094],
                    #  ["/d2/fernanda/vcf_base/c11.vcf", 4045628],
                    #  ["/d2/fernanda/vcf_base/c19.vcf", 1832506],
                    #  ["/d2/fernanda/vcf_base/c14.vcf", 2655067],
                    #  ["/d2/fernanda/vcf_base/c2.vcf", 7081600],
                    #  ["/d2/fernanda/vcf_base/c13.vcf", 2857916]
                    ]
    print(originalSizes)
    # [([id, [sizes...]), ([id, [sizes...]), ...]
    samples = 10
    newSizes = util.GenerateExpectedSizes(originalSizes, samples)
    print(newSizes)

    # for i in range(1):
    #     folder = os.path.join(destPath, str(i+1))
    #     if not os.path.isdir(folder):
    #         os.mkdir(folder)

    for i in range(len(originalSizes)):
        sizes = newSizes[i]
        name = originalSizes[i][0].split("/")[-1]
        for id, size in sizes:
            if id != 1:
                continue
            # Keep the same name, but saved in a new subfolder
            new_file = os.path.join(destPath, name)
            source = open(originalSizes[i][0], 'r')
            with open(new_file, 'w') as result:
                # # Keep metadata
                line = source.readline()
                while (line[0:2] == "##"):
                    result.write(line)
                    line = source.readline()

                # hea(der line
                result.write(line)

                mask = util.GetMask(size)
                for j in range(size):
                    if mask[j] == 1:
                        result.write(source.readline())
                    else:
                        source.readline()
            source.close()
            print(id, "finished for", name)
