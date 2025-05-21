from evaluadorBlosum import evaluadorBlosum
import copy
import math
import random

class StochasticLocalSearch:
    def __init__(self, evaluator):
        self.evaluator = evaluator
        self.max_iterations = 5  # Iteraciones máximas por bacteria
        self.max_no_improve = 2  # Iteraciones sin mejora para detenerse
        self.neighborhood_size = 3  # Tamaño del vecindario a explorar
        
    def optimize(self, bacterium, current_score):
        best_bacterium = copy.deepcopy(bacterium)
        best_score = current_score
        no_improve = 0
        
        for _ in range(self.max_iterations):
            if no_improve >= self.max_no_improve:
                break
                
            improved = False
            neighbors = self.generate_neighbors(best_bacterium)
            
            for neighbor in neighbors:
                neighbor_score = self.evaluate_bacterium(neighbor)
                if neighbor_score > best_score:
                    best_bacterium = copy.deepcopy(neighbor)
                    best_score = neighbor_score
                    improved = True
                    no_improve = 0
            
            if not improved:
                no_improve += 1
                
        return best_bacterium, best_score
        
    def generate_neighbors(self, bacterium):
        neighbors = []
        for _ in range(self.neighborhood_size):
            neighbor = copy.deepcopy(bacterium)
            # Operador 1: Mover gap aleatorio
            seq_idx = random.randint(0, len(neighbor)-1)
            if '-' in neighbor[seq_idx]:
                gap_pos = random.choice([i for i, x in enumerate(neighbor[seq_idx]) if x == '-'])
                new_pos = random.randint(0, len(neighbor[seq_idx])-1)
                neighbor[seq_idx].pop(gap_pos)
                neighbor[seq_idx].insert(new_pos, '-')
                neighbors.append(neighbor)
            
            # Operador 2: Intercambiar dos gaps de diferentes secuencias
            if len(bacterium) > 1:
                neighbor = copy.deepcopy(bacterium)
                seq1, seq2 = random.sample(range(len(neighbor)), 2)
                gaps1 = [i for i, x in enumerate(neighbor[seq1]) if x == '-']
                gaps2 = [i for i, x in enumerate(neighbor[seq2]) if x == '-']
                if gaps1 and gaps2:
                    gap1 = random.choice(gaps1)
                    gap2 = random.choice(gaps2)
                    neighbor[seq1][gap1], neighbor[seq2][gap2] = neighbor[seq2][gap2], neighbor[seq1][gap1]
                    neighbors.append(neighbor)
        
        return neighbors
        
    def evaluate_bacterium(self, bacterium):
        score = 0
        evaluator = evaluadorBlosum()
        
        for col in range(len(bacterium[0])):
            column = [bacterium[row][col] for row in range(len(bacterium))]
            pairs = self.obtener_pares_unicos(column)
            for pair in pairs:
                score += evaluator.getScore(pair[0], pair[1])
        
        return score
        
    def obtener_pares_unicos(self, columna):
        pares_unicos = set()
        for i in range(len(columna)):
            for j in range(i+1, len(columna)):
                par = tuple(sorted([columna[i], columna[j]]))
                pares_unicos.add(par)
        return list(pares_unicos)