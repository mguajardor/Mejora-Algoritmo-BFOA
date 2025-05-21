import copy
import math
from multiprocessing import Manager, Pool, managers
from pickle import FALSE, TRUE
from evaluadorBlosum import evaluadorBlosum
import numpy
from fastaReader import fastaReader
import random
from copy import copy
import copy
import concurrent.futures

class bacteria():
    def __init__(self, numBacterias):
        manager = Manager()
        self.blosumScore = manager.list(range(numBacterias))
        self.tablaAtract = manager.list(range(numBacterias))
        self.tablaRepel = manager.list(range(numBacterias))
        self.tablaInteraction = manager.list(range(numBacterias))
        self.tablaFitness = manager.list(range(numBacterias))
        self.granListaPares = manager.list(range(numBacterias))
        self.NFE = manager.list(range(numBacterias))

    def resetListas(self, numBacterias):
        manager = Manager()
        self.blosumScore = manager.list(range(numBacterias))
        self.tablaAtract = manager.list(range(numBacterias))
        self.tablaRepel = manager.list(range(numBacterias))
        self.tablaInteraction = manager.list(range(numBacterias))
        self.tablaFitness = manager.list(range(numBacterias))
        self.granListaPares = manager.list(range(numBacterias))
        self.NFE = manager.list(range(numBacterias))
        
    def cuadra(self, numSec, poblacion):
        for i in range(len(poblacion)):
            bacterTmp = poblacion[i]
            bacterTmp = list(bacterTmp)
            bacterTmp = bacterTmp[:numSec]
            maxLen = 0
            for j in range(numSec):
                if len(bacterTmp[j]) > maxLen:
                    maxLen = len(bacterTmp[j])
                    for t in range(numSec):
                        gap_count = maxLen - len(bacterTmp[t])
                        if gap_count > 0:
                            bacterTmp[t].extend(["-"] * gap_count)
                            poblacion[i] = tuple(bacterTmp)
    
    def limpiaColumnas(self):
        i = 0
        while i < len(self.matrix.seqs[0]):
            if self.gapColumn(i):
                self.deleteCulmn(i)
            else:
                i += 1
                
    def deleteCulmn(self, pos):
        for i in range(len(self.matrix.seqs)):
            self.matrix.seqs[i] = self.matrix.seqs[i][:pos] + self.matrix.seqs[i][pos+1:]

    def gapColumn(self, col):
        for i in range(len(self.matrix.seqs)):
            if self.matrix.seqs[i][col] != "-":
                return False
        return True
    
    def tumbo(self, numSec, poblacion, numGaps):
        for i in range(len(poblacion)):
            bacterTmp = poblacion[i]
            bacterTmp = list(bacterTmp)
            for j in range(numGaps):
                seqnum = random.randint(0, len(bacterTmp)-1)
                pos = random.randint(0, len(bacterTmp[seqnum]))
                part1 = bacterTmp[seqnum][:pos]
                part2 = bacterTmp[seqnum][pos:]
                temp = part1 + ["-"] + part2
                bacterTmp[seqnum] = temp
                poblacion[i] = tuple(bacterTmp)
        
    def creaGranListaPares(self, poblacion):   
        for i in range(len(poblacion)):
            pares = list()
            bacterTmp = poblacion[i]
            bacterTmp = list(bacterTmp)
            for j in range(len(bacterTmp)):
                column = self.getColumn(bacterTmp, j)
                pares = pares + self.obtener_pares_unicos(column)
            self.granListaPares[i] = pares

    def evaluaFila(self, fila, num):
        evaluador = evaluadorBlosum()
        score = 0
        for par in fila:
            score += evaluador.getScore(par[0], par[1])
        self.blosumScore[num] = score
    
    def evaluaBlosum(self):
        with Pool() as pool:
            args = [(copy.deepcopy(self.granListaPares[i]), i) for i in range(len(self.granListaPares))]
            pool.starmap(self.evaluaFila, args)

    def getColumn(self, bacterTmp, colNum):
        column = []
        for i in range(len(bacterTmp)):
            column.append(bacterTmp[i][colNum])
        return column
            
    def obtener_pares_unicos(self, columna):
        pares_unicos = set()
        for i in range(len(columna)):
            for j in range(i+1, len(columna)):
                par = tuple(sorted([columna[i], columna[j]]))
                pares_unicos.add(par)
        return list(pares_unicos)  

    def compute_diff(self, args):
        indexBacteria, otherBlosumScore, blosumScore, d, w = args
        diff = (blosumScore[indexBacteria] - otherBlosumScore) ** 2.0
        self.NFE[indexBacteria] += 1
        return d * numpy.exp(w * diff)

    def compute_cell_interaction(self, indexBacteria, d, w, atracTrue):
        with Pool() as pool:
            args = [(indexBacteria, otherBlosumScore, self.blosumScore, d, w) for otherBlosumScore in self.blosumScore]
            results = pool.map(self.compute_diff, args)
            pool.close()
            pool.join()
    
        total = sum(results)
    
        if atracTrue:
            self.tablaAtract[indexBacteria] = total
        else:
            self.tablaRepel[indexBacteria] = total
        
    def creaTablaAtract(self, poblacion, d, w):
        for indexBacteria in range(len(poblacion)):
            self.compute_cell_interaction(indexBacteria, d, w, TRUE)

    def creaTablaRepel(self, poblacion, d, w):
        for indexBacteria in range(len(poblacion)):
            self.compute_cell_interaction(indexBacteria, d, w, FALSE)
    
    def creaTablasAtractRepel(self, poblacion, dAttr, wAttr, dRepel, wRepel):
        with concurrent.futures.ThreadPoolExecutor() as executor:
            executor.submit(self.creaTablaAtract, poblacion, dAttr, wAttr)
            executor.submit(self.creaTablaRepel, poblacion, dRepel, wRepel)

    def creaTablaInteraction(self):
        for i in range(len(self.tablaAtract)):
            self.tablaInteraction[i] = self.tablaAtract[i] + self.tablaRepel[i]

    def creaTablaFitness(self):
        for i in range(len(self.tablaInteraction)):
            valorBlsm = self.blosumScore[i]
            valorInteract = self.tablaInteraction[i]
            valorFitness = valorBlsm + valorInteract
            self.tablaFitness[i] = valorFitness
    
    def getNFE(self):
        return sum(self.NFE)
        
    def obtieneBest(self, globalNFE):
        bestIdx = 0
        for i in range(len(self.tablaFitness)):
            if self.tablaFitness[i] > self.tablaFitness[bestIdx]:
                bestIdx = i
        print("-------------------   Best: ", bestIdx, " Fitness: ", self.tablaFitness[bestIdx], 
              "BlosumScore ", self.blosumScore[bestIdx], "Interaction: ", self.tablaInteraction[bestIdx], 
              "NFE: ", globalNFE)
        return bestIdx, self.tablaFitness[bestIdx]

    def replaceWorst(self, poblacion, best):
        worst = 0
        for i in range(len(self.tablaFitness)):
            if self.tablaFitness[i] < self.tablaFitness[worst]:
                worst = i
        poblacion[worst] = copy.deepcopy(poblacion[best])