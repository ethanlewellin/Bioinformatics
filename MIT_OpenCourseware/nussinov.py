import numpy as np
import pandas as pd

def write_structure(sequence, structure):
    struct = ["." for _ in range(len(sequence))]
    for s in structure:
        struct[min(s)] = "("
        struct[max(s)] = ")"
    return "".join(struct)

def traceback(i, j, structure, DP, sequence):
    #Base Case: gone through whole sequence
    if j <= i:
        return
    
    #Unpaired J case
    elif DP[i][j] == DP[i][j-1]:
        traceback(i, j-1, structure, DP, sequence)
    
    #Paired J case
    else:
        for k in [b for b in range(i, j-min_loop_length) if pairable((sequence[b], sequence[j]))]:
            if k-1 < 0:
                if DP[i][j] == DP[k+1][j-1] + 1:
                    structure.append((k,j))
                    traceback(k+1, j-1, structure, DP, sequence)
                    break
            elif DP[i][j] == DP[i][k-1] + DP[k+1][j-1] + 1:
                structure.append((k,j))
                traceback(i, k-1, structure, DP, sequence)
                traceback(k+1, j-1, structure, DP, sequence)
                break

#Check if two nucleotides can be paired
def pairable(tup):
    if tup in [('A', 'T'), ('T', 'A'), ('C', 'G'), ('G', 'C'),('A', 'U'), ('U', 'A')]:
        return True
    return False

#Function to get optimal pairs
def optimal_pairings(i,j,sequence):
    if i >= j-min_loop_length:
        return 0
    else:
        #i,j can be paired or unpaired, if unpaired recussive optimal scoring
        unpaired = optimal_pairings(i, j-1, sequence)
        
        
        pairs = [0]
        for k in range(i,j-min_loop_length):
            if pairable((sequence[j], sequence[k])):
                pairs.append(1 + optimal_pairings(i,k-1, sequence) + optimal_pairings(k+1,j-1, sequence))
        
        paired = max(pairs)
        return max(unpaired, paired)

def nussinov(sequence, min_loop_lengthVal):
    
    global min_loop_length
    min_loop_length = min_loop_lengthVal
    #Initialize NxN matrix of zeros that can't have pairings
    N = len(sequence)
    score_matrix = np.empty((N,N))
    score_matrix[:] = np.NAN
    
    for i in range (0, min_loop_length):
        for j in range(N-i):
            k = i+j
            score_matrix[j][k] = 0
            
    structure = []
    
    #fill the matrix
    for i in range(min_loop_length, N):
        for j in range(N-i):
            k = i+j
            score_matrix[j][k] = optimal_pairings(j,k,sequence)
            
    for i in range(N):
        for j in range(0,i):
            score_matrix[i][j] = score_matrix[j][i]
            
    
    traceback(0,N-1,structure,score_matrix,sequence)
    
    df = pd.DataFrame(score_matrix)
    df.columns = list(sequence)
    df.index = list(sequence)
    print(df)
    print('Structure')
    print(write_structure(sequence,structure))
    
    return