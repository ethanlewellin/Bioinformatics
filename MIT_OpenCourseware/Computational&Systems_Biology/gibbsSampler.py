#%%

# 7.91 pset 3, by Colette Picard, adapted from code by Boris Zinshteyn
import random, sys, math 
from createWebLogo import createWebLogo

# ------------------------- 
# "Sequence" class that represents sequences, with
# some fields and methods helpful for using this class in a Gibbs sampler.
# -------------------------
class Sequence:
    # fields in this class 
    # -------------------------
    seqName = ""# name of this sequence (e.g. gene name) 
    sequence = ""# nucleotide sequence 
    siteScores = []# site odds ratios =P(motif | model)/P(motif | background)
    motif = -1# current position of motif according to Gibbs sampler

    # methods in this class
    # -------------------------
    def __init__(self, name, seq):
        # initializes a new instance of a Sequence object, and initializes our
        # belief of where the motif is to be uniform across all possible motif
        # positions
        self.seqName = name 
        self.sequence = seq 
        self.siteScores = [1 for i in range(len(seq)-motifWidth+1)] + [0 for i in range(motifWidth-1)]
        self.drawNewMotifSite()

    def getMotif(self, *pos):
        # returns the motif of length motifWidth can either specify a position
        # (e.g. getMotif(4) returns the motif at position 4) or if no position
        # in specified it will return the motif at self.motif (e.g. getMotif()
        # returns self.sequence[self.motif : self.motif + motifWidth])
        if pos == (): 
            idx = self.motif 
        else: 
            idx = pos[0] 
        if idx < 0 or idx > len(self.sequence)-motifWidth: 
            print("Error - tried to access motif beginning at",idx,", index out of bounds.")
            sys.exit(1) 
        else:
            return self.sequence[idx : idx + motifWidth]

    def drawNewMotifSite(self):
        # randomly draws a new site for this Sequence's motif based on the
        # current distribution of odds ratios at each position in
        # self.siteScores INPUTS: none (uses current value of self.siteScores)
        # OUTPUTS: none (but assigns a new value of self.motif)

        # normalize the odds ratio scores to obtain a probability distribution
        tot = float(sum(self.siteScores)) 
        siteProbs = [x/tot for x in self.siteScores]

        # draw randomly according to this distribution
        randOddRatio = random.random() # returns random uniform[0,1]
        site = 0
        cumSum = siteProbs[site]
        while cumSum < randOddRatio:
            site += 1
            cumSum += siteProbs[site]
        self.motif = site

    # -------------------------

    def updateSiteScores(self, wmat, background):
        # updates the odds ratios for motifs beginning at each site in
        # self.siteScores, where odds ratio = P(motif | wmat) / P(motif |
        # background) = Pm / Pb, according to current wmat INPUTS:
        #wmat - weight matrix of current motif model in format from
        #buildWeightMatrix() background - background nucleotide frequencies
        #from findSimpleBackgroundModel()
        # OUTPUTS: none, but updates the current value of self.siteScores
        # at each site in self.sequence, calculate Pm / Pb (note that the last
        # motifWidth; 1 entries should be zero, because a motif of length
        # motifWidth cannot occur at those positions), and replace
        # self.siteScores with the new values 
        # ------------------------- 
        # 
        self.siteScores = []
        for i in range(0,len(self.sequence)):
            if i > len(self.sequence)-motifWidth:
                self.siteScores.append(0)
            else:
                Pm = 1.0 # prob under model
                Pb = 1.0 # prob under background
                for j in range(0, motifWidth):
                    Pb = Pb * background[self.getMotif(i)[j]]
                    Pm = Pm * wmat[j][self.getMotif(i)[j]]
                self.siteScores.append(Pm/Pb)
                
# ------------------------- Other helper functions -------------------------
# 

import matplotlib.pyplot as plt

def readFastaFile(filename):
    # INPUT: the name of a fasta file to read in OUTPUT: a list of Sequence
    # objects each corresponding to a sequence in the .fa file, initialized
    # according to the __init__() function under class Sequence
    try:
        fastaLines=open(filename).readlines() 
    except IOError:
        # file "filename" cannot be opened, print an error message and exit
        print("Error - could not open", filename, "\n, please check that you entered the correct file.")
        sys.exit(1) 
    
    seqList = []
    curSeq=''
    curName=None 
    
    
    for line in fastaLines:
        if '>' in line:
            curName=line.strip()[1:] 
            curSeq = ''
        else:
            curSeq=curSeq+line.strip().upper() 
            seqList.append(Sequence(curName,curSeq))
    
    return seqList


def findSimpleBackgroundModel(sequences):
    # Finds background model assuming simple 0-order model (e.g. each position
    # independent) INPUTS:
    # sequences - a list of Sequence objects
    # OUTPUTS:
    # background - a dictionary mapping the four nucleotides to their
    # frequencies across all sequences

    background = {'A': 0, 'C': 0, 'G': 0, 'T':0}
    for s in sequences:
        for nt in background:
            background[nt] = background[nt] + s.sequence.count(nt)

    # normalize
    totCounts = float(sum(background.values()))
    for nt in background:
        background[nt] = background[nt]/totCounts
    
    return background


def buildWeightMatrix(seqsToScore):
    # Builds weight matrix from motifs in all sequences except the leaveOut
    # sequence You should include pseudocounts at each position. INPUTS:
    #seqsToScore - a list of Sequence objects (left out sequence already
    #omitted)
    # OUTPUT: wmat - a list of length motifWidth, where each element corresponds
    # to a
    #position of the motif and contains a dictionary with keys = nt
    #describing the nt distribution at that position 
    # (so wmat[3]['A'] corresponds to fraction As at position 3 in motif)
    # 
    # initialize with pseudocounts at each position
    wmat = [] 
    for i in range(0, motifWidth): 
        wmat.append({'A': 1, 'C': 1, 'G':1, 'T': 1})

    # loop through all motifs, add 1 to appropriate position and nt in wmat, and
    # at the end, normalize counts to get probabilities at each position

    for s in seqsToScore:
        for i in range(0, motifWidth):
            wmat[i][s.getMotif()[i]] += 1
            
    #Normalize 
    for i in range(0, motifWidth):
        totCounts = float(sum(wmat[i].values()))
        for nt in wmat[i]:
            wmat[i][nt] = wmat[i][nt]/totCounts
    return wmat


def printWeightMatrix(wmat):
    # given a weight matrix in format given by buildWeightMatrix(), prints out
    # human-friendly version
    print("Pos\tA\tC\tG\tT" )
    for i in range(0,motifWidth): 
        print(str(i)+'\t'+str(wmat[i]['A'])+'\t'+str(wmat[i]['C'])+
              '\t'+str(wmat[i]['G'])+'\t'+str(wmat[i]['T']))


def calcRelEnt(wmat, background):
    # calculates the relative entropy of the weight matrix model wmat to the
    # background (assume every position is independent), use math.log(x,2) to
    # take the log2 of x INPUTS:
    #wmat - weight matrix in format given by buildWeightMatrix() background
    #- background nt distribution in format given by
    #findSimpleBackgroundModel()
    # OUTPUTS: the relative entropy of the weight matrix model relative to the
    # background

    relEnt = 0.0
    for pos in wmat:
        for base in pos:
            relEnt = relEnt + pos[base]*math.log(pos[base]/background[base],2)
    return relEnt


 
def getMotifScore(sequences, wmat, background):
    # the total score of a motif = sum of log2 (odds ratios for each sequence)
    score = 0 
    for s in sequences:
        # update with final weight matrix
        s.updateSiteScores(wmat, background)
        # get score at motif
        score += math.log(s.siteScores[s.motif],2) 
    return score


def printToLogo(sequences):
    # prints to the command line the motifs from each sequence in the correct
    # format for WebLogo
    seqList = []
    for s in sequences:
        seqList.append(s.getMotif())
        
    return  seqList


def gibbsSampler(fastaName, motifWidthVal, numIter):

    global motifWidth
    motifWidth = motifWidthVal
    
    #print("Running Gibbs sampler with W =", motifWidth,"for",numIter,"iterations.")

    # Read in sequences (see readFastaFile)
    sequences = readFastaFile(fastaName)
    # Find the background nucleotide distributions
    
    background = findSimpleBackgroundModel(sequences)

    # (STEP 1) pick starting sites at random done when initializing each
    # Sequence object in the readFastaFile() function - see drawNewMotifSite()
    # under Sequence class
    #print("drawNewMotifSite() not yet implemented!")
    # Repeat the following steps numIter times
    relative_entropy_list = []
    for iter in range(numIter):

        # (STEP 2) choose sequence to leave out index of sequence to be left out
        leaveOut = random.randint(0,len(sequences)-1)
        # make list of sequences with that element left out
        seqsToScore = sequences[:] 
        del seqsToScore[leaveOut]

        # (STEP 3) make weight matrix using remaining sequences
        wmat = buildWeightMatrix(seqsToScore)

        # (STEP 4) update scores across all possible motif sites for left out
        # sequence according to wmat
        sequences[leaveOut].updateSiteScores(wmat, background)

        # (STEP 5) draw a new site for motif in left out sequence at random
        # according to new distribution
        sequences[leaveOut].drawNewMotifSite()
        
        relEnt = calcRelEnt(wmat, background)
        relative_entropy_list.append(relEnt)
        
        score = getMotifScore(sequences, wmat, background)
        
    return score, wmat, sequences, background, relEnt, relative_entropy_list


def printGibbsResults(wmat, sequences, background, relEnt, relative_entropy_list):

    # Print final motif matrix, its total score and its relative entropy
    # compared to background:
    print("Final weight matrix:" )
    printWeightMatrix(wmat)
    # uncomment last line after implementing getMotifScore()
    #createWebLogo(printToLogo(sequences))
    print("Motif score =",getMotifScore(sequences, wmat, background))
    
    
    print("Final relative entropy =", relEnt)

    print('Background Distribution: \n',background)
    
    #Relitive Entropy Plot
    fig2 = plt.figure("Relative Entropy")   
    x_nums = range(1, len(relative_entropy_list) + 1) 
    plt.plot(x_nums,relative_entropy_list, 'ro') 
    plt.title('Relative Entropy of seqsA')
    plt.xlabel('Gibbs sampler iteration') 
    plt.ylabel('Relative Entropy of Motif')
    plt.show() 
# %%
