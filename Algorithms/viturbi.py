# -*- coding: utf-8 -*-
"""Viturbi.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1p-3ITt2vJ93_jTsJsjP4E6MpbJymH0_K
"""

from prompt_toolkit.layout.containers import Sequence
# pHMM

import numpy as np
import pandas as pd
import networkx.drawing.nx_pydot as gl
import networkx as nx
import matplotlib.pyplot as plt
from pprint import pprint
##matplotlib inline

#Get sequences
alignments = open('align.txt','r').read()
alignments=alignments.upper()
sequences = alignments.splitlines()
matches = sequences[-1]
sequences.pop()
sequences = list(map(list,sequences))

# Get States
states = ["B",'I0']
for i in range(matches.count('*')):
  states.append("M" + str(i+1))
  states.append("I" + str(i+1))
  states.append("D" + str(i+1))
states.append("E")

#Temp State to make looping easier
states.insert(2, "D0")

#Get transition represenation of sequences

def check_state(idx):
  if idx < matchIdx[0]:
    return 0
  else:
    for i in range(len(matchIdx)):
      if i == len(matchIdx)-1:
        break
      if idx > matchIdx[i] and idx < matchIdx[i+1]:
        return i+1
  return len(matchIdx)

matchIdx = []
for i in range(len(matches)):
    if matches[i] == "*":
      matchIdx.append(i)

tSequences = []

for j in range(len(sequences)):
  tseq = ['B']
  for i in range(len(sequences[j])):
    if not i in matchIdx and sequences[j][i] != '-':
      tseq.append("I"+str(check_state(i)))
    if i in matchIdx:
      if sequences[j][i] != '-':
        tseq.append("M"+str(matchIdx.index(i)+1))
      else:
        tseq.append("D"+str(matchIdx.index(i)+1))
  tseq.append('E')
  tSequences.append(tseq)

#display(tSequences)

#Get Transition Probabilities
trMat = np.zeros((len(states),len(states)))
for i in range(len(tSequences)):
  for j in range(len(tSequences[i])):
    if tSequences[i][j] != "E":
      trMat[states.index(tSequences[i][j]),states.index(tSequences[i][j+1])] += 1

for i in range(len(trMat)):
  if sum(trMat[i]) != 0:
    trMat[i] = trMat[i]/sum(trMat[i])

display(pd.DataFrame(trMat, columns = states, index=states))

# Get match stat probabilities
dna = ["A","T","C","G"]

msProb = np.zeros((len(matchIdx),len(dna)))
for i in range(len(matchIdx)):
  for j in range(len(sequences)):
    if sequences[j][matchIdx[i]] in dna:
      msProb[i,dna.index(sequences[j][matchIdx[i]])] += 1

for i in range(len(msProb)):
  if sum(msProb[i]) != 0:
    msProb[i] = msProb[i]/sum(msProb[i])

display(pd.DataFrame(msProb, columns = dna))

from numpy.core.fromnumeric import shape
import graphviz
G = graphviz.Digraph('morkov',engine="neato")
G.graph_attr["rankdir"] = "LR"

# nodes correspond to states
Ms = []
Is = ['I0']
Ds = []
for i in range(len(matchIdx)):
  Ms.append("M"+str(i+1))
  Is.append("I"+str(i+1))
  Ds.append("D"+str(i+1))

#Add begin node
G.node(name = "B", label = "B", pos='0,0!')

#Add insertion nodes
dist=0
for i in range(len(Is)):
  G.node(Is[i], label = Is[i], pos=str(dist+2.5)+',1.5!')
  dist+=2.5

#Add deletion nodes
dist=0
for i in range(len(Ds)):
  if Ds[i] == 'D0':
    continue #Don't plot temp node
  G.node(Ds[i], label = Ds[i], pos = str(dist+2.5)+',2.75!')
  dist+=2.5

#Add match nodes
dist=0
for i in range(len(Ms)):
  G.node(Ms[i], label = "{" + "M"+ str(i+1) + "|{<p1> " + dna[0] + ": " + str(round(msProb[i,0],3)) + "|<p2>     " + dna[1] + ": " + str(round(msProb[i,1],3)) +"|<p3>     " + dna[2] + ": " + str(round(msProb[i,2],3)) + "|<p4>     " + dna[3] + ": " + str(round(msProb[i,3],3)) +"} }",
         shape = "record",pos = str(dist+2.5)+',0!')
  dist+=2.5

#Add end node
G.node(name = "E", label = "E",pos=str(dist+2.5)+',0!')
#print(G.source)

#get edges from trasition probabillity matrix

def get_markov_edges(transition_Matrix):
    trDf = pd.DataFrame(transition_Matrix, columns=states, index=states)
    edges = {}
    for col in trDf.columns:
        for idx in trDf.index:
            edges[(idx,col)] = trDf.loc[idx,col]
    return edges

edges_wts = get_markov_edges(trMat)
dic_out = {}
for x, y in edges_wts.items():
    if y != 0:
        dic_out[x] = round(y,3)

edges_wts = dic_out
#pprint(edges_wts)

for k, v in edges_wts.items():
    tmp_origin, tmp_destination = k[0], k[1]
    G.edge(k[0],k[1],label= str(v))

#print(G.source)

#Render picture
G.render('markov.gv').replace('\\', '/')
G.render('doctest-output/markov.gv', view=True, format = "jpg").replace('\\', '/')

from pandas.core.arrays.sparse import array
#Viturbi Algorithm

Seq = 'AATGAC'
Seq = '-'+Seq
Seq = list(Seq)

viturbi = np.zeros((len(states), len(Seq)))
dirMat = np.repeat('--',len(states)*len(Seq)).reshape(len(states),len(Seq))

#Set begin probability to 1
viturbi[0,0] = 1

#Start Filling Matricies
for j in range(len(Seq)):
  for i in range(len(states)):

    if j>0:
      # Match States
      if i%3 == 0 and i>0:
        #from previous match
        M = viturbi[i-3,j-1]*trMat[i-3,i]
        #from previous insert
        I = viturbi[i-2,j-1]*trMat[j-2,i]
        #from previous deletion
        D = viturbi[i-1,j-1]*trMat[j-2,i]

        #print(M,I,D)
        max_value = max(M,I,D)
        if(states[i] != 'E'):
          #print(states[i])
          #print(i,int((i-3)/3))
          viturbi[i,j] = max_value*msProb[int((i-3)/3),dna.index(Seq[j])]
        else:
          viturbi[i,j] = max(M,I)

        if M == max_value:
          dirMat[i,j] = states[i-3]
        elif I == max_value:
          dirMat[i,j] = states[i-2]
        else:
          dirMat[i,j] = states[i-1]

        # Insertion States
      if i%3 == 1:
        #from previous match
        M = viturbi[i-1,j-1]*trMat[i-1,i]
        #from previous delete
        D = viturbi[i-2,j-1]*trMat[i-2,i]
        #from self
        I = viturbi[i,j-1]*trMat[i,i]

        #print(M,I,D)
        max_value = max(M,I,D)
        viturbi[i,j] = max_value*0.25

        if M == max_value:
          dirMat[i,j] = states[i-1]
        elif I == max_value:
          dirMat[i,j] = states[i]
        else:
          dirMat[i,j] = states[i-2]


    # Delete States
    if i%3 == 2:
      #from previous match
      M = viturbi[i-5,j]*trMat[i-5,i]
      #from previous insert
      I = viturbi[i-4,j]*trMat[i-4,i]
      #from previous delete
      D = viturbi[i-3,j]*trMat[i-3,i]

      #print(i,j,max(M,I,D))
      max_value = max(M,I,D)
      viturbi[i,j] = max_value


      if M == max_value:
        dirMat[i,j] = states[i-5]
      elif I == max_value:
        dirMat[i,j] = states[i-4]
      else:
        dirMat[i,j] = states[i-3]

states.remove('D0')
viturbi = np.delete(viturbi,2,axis=0)
dirMat = np.delete(dirMat,2, axis=0)

display(pd.DataFrame(viturbi, columns=Seq, index=states))
#pd.DataFrame(dirMat, columns=Seq, index=states)

from tables import index
from numpy.lib.function_base import append
# Trace Back

i=len(states)-1
j=len(Seq)-1
MLStates = ['E']
while(dirMat[i,j]!='--'):
  MLStates.append(dirMat[i,j])
  i = states.index(dirMat[i,j])
  if list(dirMat[i,j])[0] != 'D':
    j = j-1

MLStates.reverse()
MLStates

