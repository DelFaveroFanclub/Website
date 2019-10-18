
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 12 14:00:48 2019

@author: Michiel
"""

from .Gene import BED, Motif

#Function accepting webtool input
def FindMotif(file = '',motif=Motif(), overlap = 0.8):
    bed = BED(file)
    genelist = bed.get_sequences()
    frequent_patterns = []
    for gene in genelist:
         presence = motif.is_present_in(gene,min_overlap=overlap)
         frequent_patterns.append(presence)
    frequency = sum(frequent_patterns)/len(frequent_patterns)

    motif_len = len(motif)
    motif_num = len(motif.all_motifs())
    total_len = sum([len(gene) for gene in genelist]) - motif_len*len(genelist)
    exp_freq = 1- (1-1/4**motif_len)**(total_len*motif_num)
    return (frequency,exp_freq)

#Helper frequent pattern miner, currently not used
def findfrequentpattern(string):
    charcounter = dict()

    charlist = list(filter(lambda x: x in 'AGTC', string))
    for char in charlist:
        if char in charcounter:
            charcounter[char] += 1
        else:
            charcounter[char] = 1

    maxchar = max(charcounter)

    return (maxchar, charcounter[char])


#For testing
if __name__ == '__main__':
    print(FindMotif('bed_ex'))
