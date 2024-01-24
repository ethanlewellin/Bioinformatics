#%%
from collections import Counter
import math
import matplotlib.pyplot as plt 

def createWebLogo(seqList):    

    # Create a sequence logo plot with cumulative entropy for each nucleotide
    positions = range(1, len(seqList[0]) + 1)  # Assuming all sequences have the same length

    total_bases = len(seqList)
    # Calculate the maximum possible entropy
    max_possible_entropy = math.log2(len('ACGT'))
    
    # Initialize lists to store data for bar chart

    # Create a bar chart
    fig, ax = plt.subplots()
    handles = []  # To store legend handles
    #Define colors for each nucleotide
    nucleotide_colors = {'A': 'green', 'C': 'blue', 'G': 'orange', 'T': 'red'}

    # Iterate over each position
    for pos in positions:
        # Calculate base frequencies at the current position
        base_counts = Counter(seq[pos - 1] for seq in seqList)
        
        #print(base_counts)
        # Initialize a dictionary to store cumulative entropy for each nucleotide
        difference_entropy = {base: 0 for base in 'ACGT'}
        observed_probability = {base: 0 for base in 'ACGT'}
        # Calculate the observed entropy for each nucleotide
        for base, count in base_counts.items():
            observed_probability[base] += count / total_bases
            entropy_term = -(observed_probability[base] * math.log2(observed_probability[base])) if observed_probability[base] > 0 else 0
            difference_entropy[base] += entropy_term
        errorCorrectionFactor = (len('ATCG') - 1) / (2 * math.log2(2) * len(seqList))

        #Get heights for graphing
        base_heights = {base: 0 for base in 'ACGT'}
        for base in observed_probability.keys():
            height = observed_probability[base] * (max_possible_entropy - (sum(difference_entropy.values()) + errorCorrectionFactor))
            base_heights[base] += height if height > 0 else 0
        # Order the nucleotides based on cumulative entropy from least to most
        ordered_bases = sorted(base_heights.keys(), key=lambda x: base_heights[x])
                
        # Store data for bar chart
        #for i in positions:
        bottom = 0
        if pos == 1:
            handle_labels = ordered_bases
        for base in ordered_bases:
            color = nucleotide_colors[base]  # Use the color corresponding to the nucleotide
            bar = ax.bar(pos, base_heights[base], label=base, bottom=bottom, color=color)
            if base_heights[base] > max_possible_entropy/ (8*3):
                ax.bar_label(bar, label_type='center', labels = base)
            handles.append(bar[0])  # Store the first element of each bar as a handle
            bottom += base_heights[base]

    plt.xlabel('Position')
    plt.ylabel('Bits')
    plt.ylim(0, max_possible_entropy)
    plt.title('Sequence Logo with Bits')
    plt.xticks(positions)
    plt.legend(handles=handles, labels=handle_labels)
    
    return

createWebLogo(['ATCG','TTAG','CACC'])
# %%
