from ast import If
import sys
import os
import numpy as np
from random import randint, random
from VCFUtils import VCFUtils

if __name__ == "__main__":
    reference = sys.argv[1]
    destPath = sys.argv[2]
    util = VCFUtils()
    originalNames = []

    # lens = [4, 6, 8, 10, 12, 14, 16, 18, 20]
    lens = [4, 6, 8]
    samples = 100

    with open(reference, 'r') as VCFList:
        # Recuperamos las lineas de la referencia (solo 1 cromosoma)
        lines = VCFList.readlines()[1:]

    # Por cada largo de patron deseado
    for l in lens:
        with open(os.path.join(destPath, "p" + str(l) + ".txt"), 'w') as result:
            print("Inicia recoleccion de patrones de largo", l)
            # Para todas las lineas posibles generamos un arreglo aleatorio para ver si usamos la linea o no
            doIt = np.zeros(len(lines))
            doIt[:samples * 2] = 1
            done = 0
            np.random.shuffle(doIt)
            # Por cada linea
            for i in range(len(lines)):
                if done == samples: break
                # Si se usa
                if doIt[i]:
                    line = lines[i]
                    start = randint(0, len(line) - l - 1)
                    if line[start:start+l] == "A"*l: continue
                    done+=1
                    # print(line[start:start+l], i)
                    result.write(line[start:start+l] + "\n")
            print(" Finaliza recoleccion de patrones de largo", l)
            
            

