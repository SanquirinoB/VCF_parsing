from concurrent.futures import process
import sys
import os
import re
from VCFUtils import VCFUtils

if __name__ == "__main__":
    processType = sys.argv[1]
    source = sys.argv[2]
    destPath = sys.argv[3]
    util = VCFUtils()
    
    # Generamos una referencia nueva en el folder de destino
    if (processType == "REF"):
        aux = source.split("/")[-1]
        new_reference = os.path.join(destPath, "ref.fa")
        refFile = open(source, 'r')
        print("Start", source , "->", new_reference)
        with open(new_reference, 'w') as newFile:
            line = refFile.readline()
            while(line):
                if(line[0] != ">"):
                    line = re.sub('[^ACGT\n]', 'A', line)
                newFile.write(line)
                line = refFile.readline()

        refFile.close()
    # Generamos un VCF limpio y eliminamos el anterior
    elif (processType == "VCF"):
        originalNames = []
        # Recuperamos la lista de nombres de los archivos a prcoesar
        with open(source, 'r') as VCFList:
            line = VCFList.readline()
            print(line)
            while line:
                originalNames.append(line[:-1])
                line = VCFList.readline()

        for name in originalNames:
            aux = name.split("/")[-1]
            newName = os.path.join(destPath, "s100", aux)
            print("Start", name , "->", newName)
            VCF = open(name, 'r')
            with open(newName, 'w') as newVCF:

                # Mantenemos la meta data
                line = VCF.readline()
                while (line[0:2] == "##"):
                    newVCF.write(line)
                    line = VCF.readline()

                # Mantenemos el header
                newVCF.write(line)
                # Por cada registro eliminaremos forzosamente la letra N
                line = VCF.readline()

                while(line):
                    if (len(line) == 2): break
                    line = line.split("\t")
                    # REF
                    line[3] = re.sub('[^ACGT]', 'A', line[3])
                    # ALT
                    line[4] = re.sub('[^ACGT]', 'A', line[4])

                    newVCF.write("\t".join(line))

                    line = VCF.readline()

            print( newName,"generated")
            VCF.close()
