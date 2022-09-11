import sys
import os
from VCFUtils import VCFUtils

if __name__ == "__main__":
    vcfList = sys.argv[1]
    destPath = sys.argv[2]
    util = VCFUtils()
    originalNames = []
    # Recuperamos la lista de nombres de los archivos a prcoesar
    with open(vcfList, 'r') as VCFList:
        line = VCFList.readline()
        print(line)
        while line:
            originalNames.append(line[:-1])
            line = VCFList.readline()

    # Ver que todo esta en orden
    print(originalNames)

    # Recuperar cantidad de samples en estos VCF a partir solo del primer archivo
    originalSize = util.CountSamples(originalNames[0])
    # originalSizes = util.CountSamples(originalNames[0])

    # Ver que todo esta en orden
    print(originalSize)

        # Esta linea es para generar el primer sampling sobre los VCF original
    samples = 1
        # Despues usamos esta para escalar sobre el 10% original
    # samples = 10
    newSizes = util.GenerateExpectedSizes(originalSize, samples)
    print(newSizes)

    # Creo la carpeta donde almacenare los archivos generados
    for i in range(samples):
        folder = os.path.join(destPath, "s" + str((i+1) * 10))
        if not os.path.isdir(folder):
            os.mkdir(folder)

    # Por cada archivo
    for i in range(len(originalNames)):
        name = originalNames[i].split("/")[-1]
        # Por cada tamano a generar
        for id, size in newSizes:
            # Creamos el nuevo documento, mismo nombre pero diferente path
            new_file = os.path.join(destPath, "s" + str((id) * 10), name)
            # Abrimos el source
            source = open(originalNames[i], 'r')
            # Empezamos a escribir el nuevo documento
            with open(new_file, 'w') as result:
                # Mantenemos la meta data
                line = source.readline()
                while (line[0:2] == "##"):
                    result.write(line)
                    line = source.readline()

                # Desde la header line, empezaremos a conservar solo los samples esperados
                aux = line.split("\t")[:9 + size]
                aux.append("\n")
                result.write("\t".join(aux))
                line = source.readline()

                while(line):
                    aux = line.split("\t")[:9 + size]
                    aux.append("\n")
                    result.write("\t".join(aux))
                    line = source.readline()

                aux = line.split("\t")[:9 + size]
                aux.append("\n")
                result.write("\t".join(aux))

            source.close()
            print(id, "finished for", name)
