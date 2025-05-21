from bacteria import bacteria
from fastaReader import fastaReader
from evaluadorBlosum import evaluadorBlosum
from StochasticLocalSearch import StochasticLocalSearch
import time
import copy
from multiprocessing import Manager
import random

if __name__ == "__main__":
    numeroDeBacterias = 4
    numRandomBacteria = 1
    iteraciones = 3
    tumbo = 2
    nado = 3
    secuencias = list()
    dAttr = 0.1
    wAttr = 0.002
    hRep = dAttr
    wRep = .001
    
    secuencias = fastaReader().seqs
    names = fastaReader().names
    
    for i in range(len(secuencias)):
        secuencias[i] = list(secuencias[i])

    globalNFE = 0

    
    manager = Manager()
    numSec = len(secuencias)
    print("numSec: ", numSec)
    
    poblacion = manager.list(range(numeroDeBacterias))
    names = manager.list(names)
    NFE = manager.list(range(numeroDeBacterias))

    def poblacionInicial():
        for i in range(numeroDeBacterias):
            bacterium = []
            for j in range(numSec):
                bacterium.append(secuencias[j])
            poblacion[i] = list(bacterium)

    def printPoblacion():
        for i in range(numeroDeBacterias):
            print(poblacion[i])
            
    operadorBacterial = bacteria(numeroDeBacterias)    
    veryBest = [None, None, None]

    start_time = time.time()
    
    print("poblacion inicial ...")
    poblacionInicial() 
    
    for it in range(iteraciones):
        print("poblacion inicial creada - Tumbo ...")
        operadorBacterial.tumbo(numSec, poblacion, tumbo)
        print("Tumbo Realizado - Cuadrando ...")
        operadorBacterial.cuadra(numSec, poblacion)
        print("poblacion inicial cuadrada - Creando granLista de Pares...")
        operadorBacterial.creaGranListaPares(poblacion)
        print("granList: creada - Evaluando Blosum Parallel")
        operadorBacterial.evaluaBlosum()
        print("blosum evaluado - creando Tablas Atract Parallel...")
        operadorBacterial.creaTablasAtractRepel(poblacion, dAttr, wAttr, hRep, wRep)
        operadorBacterial.creaTablaInteraction()
        print("tabla Interaction creada - creando tabla Fitness")
        operadorBacterial.creaTablaFitness()
        print("tabla Fitness creada")
        globalNFE += operadorBacterial.getNFE()
        
        # NUEVA PARTE: Búsqueda Local Estocástica
        print("Aplicando búsqueda local estocástica...")
        sls = StochasticLocalSearch(evaluadorBlosum())
        
        for i in range(len(poblacion)):
            current_bacterium = copy.deepcopy(poblacion[i])
            current_score = operadorBacterial.tablaFitness[i]
            
            improved_bacterium, improved_score = sls.optimize(current_bacterium, current_score)
            
            if improved_score > current_score:
                poblacion[i] = improved_bacterium
                operadorBacterial.tablaFitness[i] = improved_score
                operadorBacterial.blosumScore[i] = sls.evaluate_bacterium(improved_bacterium)
                globalNFE += sls.max_iterations * sls.neighborhood_size
        
        bestIdx, bestFitness = operadorBacterial.obtieneBest(globalNFE)
        if (veryBest[0] is None) or (bestFitness > veryBest[1]):
            veryBest[0] = bestIdx
            veryBest[1] = bestFitness
            veryBest[2] = copy.deepcopy(poblacion[bestIdx])
        operadorBacterial.replaceWorst(poblacion, veryBest[0])
        operadorBacterial.resetListas(numeroDeBacterias)

    print("Very Best: ", veryBest)
    print("--- %s seconds ---" % (time.time() - start_time))