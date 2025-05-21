import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import time
from datetime import datetime
import os
import copy
from multiprocessing import Manager
from bacteria import bacteria
from fastaReader import fastaReader
from evaluadorBlosum import evaluadorBlosum
from StochasticLocalSearch import StochasticLocalSearch

def ejecutar_algoritmo_mejorado():
    """Ejecuta una corrida del algoritmo mejorado con búsqueda local"""
    # Configuración de parámetros
    numeroDeBacterias = 6
    iteraciones = 3
    tumbo = 2
    dAttr = 0.1
    wAttr = 0.002
    hRep = 0.1
    wRep = 0.001
    
    # Inicialización
    secuencias = [list(seq) for seq in fastaReader().seqs]
    manager = Manager()
    poblacion = manager.list(range(numeroDeBacterias))
    operador = bacteria(numeroDeBacterias)
    sls = StochasticLocalSearch(evaluadorBlosum())
    veryBest = [None, None, None]
    globalNFE = 0
    
    # Población inicial
    for i in range(numeroDeBacterias):
        poblacion[i] = copy.deepcopy(secuencias)
    
    # Datos para registro
    fitness_history = []
    tiempo_inicio = time.time()
    
    for it in range(iteraciones):
        # Fase normal del algoritmo
        operador.tumbo(len(secuencias), poblacion, tumbo)
        operador.cuadra(len(secuencias), poblacion)
        operador.creaGranListaPares(poblacion)
        operador.evaluaBlosum()
        operador.creaTablasAtractRepel(poblacion, dAttr, wAttr, hRep, wRep)
        operador.creaTablaInteraction()
        operador.creaTablaFitness()
        globalNFE += operador.getNFE()
        
        # Fase de búsqueda local mejorada
        for i in range(len(poblacion)):
            current = copy.deepcopy(poblacion[i])
            current_score = operador.tablaFitness[i]
            improved, improved_score = sls.optimize(current, current_score)
            if improved_score > current_score:
                poblacion[i] = improved
                operador.tablaFitness[i] = improved_score
                operador.blosumScore[i] = sls.evaluate_bacterium(improved)
                globalNFE += sls.max_iterations * sls.neighborhood_size
        
        # Registro de resultados
        best_idx, best_fitness = operador.obtieneBest(globalNFE)
        fitness_history.append(best_fitness)
        
        if veryBest[0] is None or best_fitness > veryBest[1]:
            veryBest = [best_idx, best_fitness, copy.deepcopy(poblacion[best_idx])]
            
        operador.replaceWorst(poblacion, veryBest[0])
        operador.resetListas(numeroDeBacterias)
    
    tiempo_ejecucion = time.time() - tiempo_inicio
    return {
        'fitness_history': fitness_history,
        'best_fitness': veryBest[1],
        'execution_time': tiempo_ejecucion,
        'nfe': globalNFE,
        'best_alignment': veryBest[2]
    }

def ejecutar_30_corridas_mejorado():
    """Ejecuta 30 corridas del algoritmo mejorado y genera gráficas"""
    # Crear directorio para resultados si no existe
    if not os.path.exists('resultados_mejorado'):
        os.makedirs('resultados_mejorado')
    
    resultados = []
    for i in range(30):
        print(f"Ejecutando corrida {i+1}/30 del algoritmo mejorado...")
        resultados.append(ejecutar_algoritmo_mejorado())
    
    # Crear dataframe con resultados
    df = pd.DataFrame(resultados)
    
    # Guardar resultados en CSV
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    df.to_csv(f'resultados_mejorado/resultados_{timestamp}.csv', index=False)
    
    # Guardar el mejor alineamiento encontrado
    mejor_corrida = df['best_fitness'].idxmax()
    mejor_alineamiento = resultados[mejor_corrida]['best_alignment']
    with open(f'resultados_mejorado/mejor_alineamiento_{timestamp}.fasta', 'w') as f:
        for i, seq in enumerate(mejor_alineamiento):
            f.write(f">Secuencia_{i+1}\n")
            f.write(''.join(seq) + '\n')
    
    # Generar gráficas
    generar_graficas_mejorado(df, timestamp)
    
    return df

def generar_graficas_mejorado(df, timestamp):
    """Genera gráficas para el algoritmo mejorado"""
    plt.figure(figsize=(15, 10))
    
    # Gráfica de convergencia de todas las corridas
    plt.subplot(2, 2, 1)
    for i, history in enumerate(df['fitness_history']):
        plt.plot(history, alpha=0.5, label=f'Corrida {i+1}' if i < 5 else None)
    plt.title('Convergencia de todas las corridas')
    plt.xlabel('Iteración')
    plt.ylabel('Fitness')
    if len(df) <= 5:
        plt.legend()
    plt.grid()
    
    # Gráfica de convergencia promedio
    plt.subplot(2, 2, 2)
    max_len = max(len(x) for x in df['fitness_history'])
    padded = [np.pad(x, (0, max_len - len(x)), mode='edge') for x in df['fitness_history']]
    promedio = np.mean(padded, axis=0)
    std = np.std(padded, axis=0)
    plt.plot(promedio, label='Promedio')
    plt.fill_between(range(len(promedio)), promedio-std, promedio+std, alpha=0.2)
    plt.title('Convergencia Promedio ± Desviación Estándar')
    plt.xlabel('Iteración')
    plt.ylabel('Fitness')
    plt.grid()
    
    # Gráfica de distribución del mejor fitness
    plt.subplot(2, 2, 3)
    plt.boxplot(df['best_fitness'])
    plt.title('Distribución del Mejor Fitness')
    plt.ylabel('Fitness')
    plt.grid()
    
    # Gráfica de tiempo de ejecución vs fitness
    plt.subplot(2, 2, 4)
    plt.scatter(df['execution_time'], df['best_fitness'])
    plt.title('Tiempo de Ejecución vs Mejor Fitness')
    plt.xlabel('Tiempo (s)')
    plt.ylabel('Fitness')
    plt.grid()
    
    plt.tight_layout()
    plt.savefig(f'resultados_mejorado/graficas_{timestamp}.png')
    plt.show()
    
    # Gráfica adicional de NFE vs fitness
    plt.figure(figsize=(8, 5))
    plt.scatter(df['nfe'], df['best_fitness'])
    plt.title('Evaluaciones de Función (NFE) vs Mejor Fitness')
    plt.xlabel('NFE')
    plt.ylabel('Fitness')
    plt.grid()
    plt.savefig(f'resultados_mejorado/nfe_vs_fitness_{timestamp}.png')
    plt.show()

if __name__ == "__main__":
    # Ejecutar las 30 corridas del algoritmo mejorado
    resultados_df = ejecutar_30_corridas_mejorado()
    
    # Mostrar resumen estadístico
    print("\nResumen estadístico del algoritmo mejorado:")
    print(resultados_df[['best_fitness', 'execution_time', 'nfe']].describe())