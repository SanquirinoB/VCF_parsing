from concurrent.futures import process
import sys
import os
from VCFUtils import VCFUtils

if __name__ == "__main__":
    processType = sys.argv[1]
    source = sys.argv[2]
    destPath = sys.argv[3]
    util = VCFUtils()
    
    # Generamos una referencia nueva en el folder de destino
    if (processType == "R"):
        new_reference = os.path.join(destPath, source)
        refFile = open(source, 'r')
        with open(new_reference, 'w') as newFile:
            line = refFile.readline()
            while(line):
                line = line.replace("N", "A")
                line = line.replace("n", "A")
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
            aux = name.split(".")
            newName = "".join(aux[0]) + "_noN.vcf"

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
                line = VCF.readline().split("\t")
                while(line):
                    # REF
                    line[3] = line[3].replace("N", "A")
                    # ALT
                    line[4] = line[4].replace("N", "A")
                    newVCF.write("\t".join(line))
                    line = VCF.readline().split("\t")
            VCF.close()
            os.remove(name)
