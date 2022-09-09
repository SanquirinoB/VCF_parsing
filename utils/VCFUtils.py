import numpy as np

class VCFUtils():

    def __init__(self) -> None:
        pass

    def GetSizes(self, nameLists):
        n_record = 0
        sizePerFile = []
        for name in nameLists:
            with open(name, 'r') as VCF:
                line = VCF.readline()
                
                # Skip metadata
                while (line[0:2] == "##"):
                    line = VCF.readline()

                while (line):
                    n_record += 1
                    line = VCF.readline()

                print(name, n_record)
                sizePerFile.append(n_record)
            n_record = 0

        return sizePerFile

    def GenerateExpectedSizes(self, originalSizes, n = 10):

        newSizes = []
        for _, value in originalSizes:
            aux = []
            for i in range(n):
                aux.append([ i + 1, int(value * ((i + 1) / 10))])
            newSizes.append(aux)

        return newSizes

    def GetMask(self, n):
        return list(np.random.randint(2, size=(n)))
