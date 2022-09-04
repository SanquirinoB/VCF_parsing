import sys
from VCFUtils import VCFUtils

if __name__ == "__main__":
    vcfList = sys.argv[1]
    destPath = sys.argv[2]
    util = VCFUtils()
    originalNames = []

    with open(vcfList, 'r') as VCFList:
        originalNames.append(VCFList.readline())

    print(originalNames)
    originalSizes = util.GetSizes(originalNames)
    print(originalSizes)
    # [([id, [sizes...]), ([id, [sizes...]), ...]
    newSizes = util.GenerateExpectedSizes(originalSizes)
    print(newSizes)
