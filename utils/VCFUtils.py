import numpy as np

class VCFUtils():

    def __init__(self) -> None:
        pass

    def CountSamples(self, vcfPath):        
        with open(vcfPath, 'r') as VCF:
            line = VCF.readline()

            while line[0:2] == "##":
                line = VCF.readline()
            
            header = line.split("\t")[9:]

        return len(header)

    def GenerateExpectedSizes(self, originalSize, n = 10):
        # Genera un escalamiento al [0.1, 0.2, ..., 0.9, 1]
        newSizes = []
        for i in range(n):
            newSizes.append([ i + 1, int(originalSize * ((i + 1) / 10))])

        return newSizes
