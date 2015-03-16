#!/usr/bin/python

"""Final Project for introscicomp2014: Antibody Isotype calculator.  Problem: our lab uses a very slow and inefficient program to determine antibody sequences.  Other labmates had developed means to analyze mutation and gene repertoire information but currently there is no way to properly determine the antibody isotype.  This program is designed to scan FASTA sequence files of immunologobulin genes and provide a best match, by percentage, to the antibody isotype that is likely being utilized. IgM, IgG1, IgG2, IgG3, IgG4, IgA1, IgA2, IgD, and IgE are all represented in our readout"""

__author__ = 'Noel Pauli (npauli@uchicago.edu)'
__version__ = '0.1.0'

#imports
import sys
import os
import glob

# this program will utilize Numpy and Biopython to do the intial sequence analysis

# constants

# functions



     # These are the two sequences to match


motifs = {'IgG1': ['CCTCCACCAAGGGCCCATCGGTCTTCCCCCTGGCACCCTCCT', 0], 
          'IgG2': ['CCTCCACCAAGGGCCCATCGGTCTTCCCCCTGGCGCCCTGCT', 0], 
          'IgG3': ['CTTCCACCAAGGGCCCATCGGTCTTCCCCCTGGCGCCCTGCT', 0], 
          'IgG4': ['CTTCCACCAAGGGCCCATCCGTCTTCCCCCTGGCGCCCTGCT', 0], 
          'IgA1': ['CATCCCCGACCAGCCCCAAGGTCTTCCCGCTGAGCCTCTGCA', 0], 
          'IgM': ['GGAGTGCATCCGCCCCAACCCTTTTCCCCCTCGTCTCCTGTG', 0], 
          'IgA2': ['CATCCCCGACCAGCCCCAAGGTCTTCCCGCTGAGCCTCGACA', 0], 
          'IgD': ['CACCCACCAAGGCTCCGGATGTGTTCCCCATCATATCAGGGT', 0], 
          'IgE': ['CCTCCACACAGAGCCCATCCGTCTTCCCCTTGACCCGCTGCT', 0]} 
   


# assign the longest sequence
# to s1, and the shortest to S2
# L1 is the length of the longest,
# l2 that of the shortest



#def sequence_length(seq1):
#    l1 = len(seq1)
#    l2 = len(motifs.get("IgG1")[0])
    #len simply counts the number of characters in the string
    #print motifs.get('IgG1')

#    if l1 >= l2:
#        s1 = seq1
#        s2 = str(motifs)
#    return indicate_startpoint(l1, l2, s1, s2)
    

# If the length of seq1 is greater than seq2, do basically nothing

# we need to find a place to start the sequence alignment. This function allows the code to iterate through the number of sequence positions on the fasta file incrimenatally.  It will stop always 42 base pairs from the end to eliminate defintion errors. 

#def indicate_startpoint(l1, l2, s1, s2):
#    for q in range(l1-l2):
#        #startpoint = z
#        score = calculate_match(s1, s2, l1, l2, q)
#        print score
#    return score
    
# function that computes a score
# by returning the number of matches
# starting from arbitrary startpoint

#finds best alignment for an Ig sequence
def calculate_match(s1, s2, l1, l2, startpoint):
    # startpoint is the point at which we want to start
    # This section is doing the actual aligning of the motif with the sequence
    matched = "" # contains string from alignement
    score = 0
    nonmatched = 0
    for i in range(l2):
        if (i + startpoint) < l1:
                # if it's matching the character
            if s1[i + startpoint] == s2[i]:
                matched = matched + "*"
                score = score + 1
            else:
                matched = matched + "-"
                nonmatched = nonmatched + 1
        
    # build some formatted output
                firstline = "." * startpoint + matched
               # print "." * startpoint + s2
               # print s1
               # print score
               # print ""
    # this is despensible once the program works, it clutters the output
  

    return (score, float(score)/(score + nonmatched))

#calculates the best score
def calculate_score(s1, motif):

    l1 = len(s1)
    l2 = len(motif)

# now try to find the best match (highest score).  Arbitrary initial score assignments
    my_best_align = None
    my_best_score = -1
    my_best_ratio = 0

#This looks for the best match in the particular motif from the sequence length minus 42bp
    for i in range(l1 - l2):
        z = calculate_match(s1, motif, l1, l2, i)
        #print z
        if z [0] > my_best_score:
            my_best_align = "." * i + motif
            my_best_score = z [0]
            my_best_ratio = z [1]
            # z is a tuple in this case with two parameters. 
        elif z [0] == my_best_score: 
            if z [1] > my_best_ratio:
                my_best_align = "." * i + motif
                my_best_score = z [0]
                my_best_ratio = z [1]   

        #print my_best_align
        #print s1
        #print "Best score:", my_best_score
       # print "Best Ratio:", my_best_ratio

#returns the best score, align, and ratio to the higher functions
    return (my_best_score, my_best_ratio, my_best_align)

def find_best_motif(s1):
# This function will iterate through the motifs to find the best match of all matches 
    my_best_motif = None
    my_best_ratio = 0
    my_best_score = -1

    
    for motif in motifs:
# pass calculate score a sequence not a name (IgG).  Get the sequence from x
        Ig_seq = motifs.get(motif)[0]
        y = calculate_score(s1, Ig_seq)

        if y[0] > my_best_score:
            print "INSIDE"
            my_best_score = y[0]
            my_best_ratio = y[1]
            my_best_motif = motif
        elif y[0] == my_best_score:
            if y[1] > my_best_ratio:
                my_best_score = [0]
                my_best_ratio = [1]
                my_best_motif = [2]

# prints the absolute winners best score, ratio, and isotype for that sequence
        print "***Best Motif: ", my_best_motif
        print "Best Score: ", my_best_score
        print "Best Ratio: ", my_best_ratio

    return (my_best_motif, my_best_score, my_best_ratio)

# import FASTA
# call find best motif
# import CSV

if __name__ == "__main__":
    seq1 = "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNATNCNGGTAGCTGGGGGGTCCTCAGGCTCTCCTGTGTAGTCTCTGGACTCACTTTCGACAATGCCTGGATGAGCTGGGTCCGCCAGGCTCCAGGGAAGGGGCTGGAGTGGGTTGGCCGTATTAAAAGCAAATCTGATGGTGGGACAACAGACTTCGCTGCACCCGTGAAAGGCAGATTCACCATCTCTAGAGATGACTCAAAAAATACAGTGTTTCTGCAAATGAACAGCCTGCAGACCGAGGACACAGCCGTGTATTACTGTGCCACAGCCCCCGGATTCCACTATTATGCTCCCTTTGACTACTGGGGCCCGGGAACCCTGGTCACCGTCTCCTCAGCCTCCACCAAGGGCCCATCGGTCTTCCCCCTGGCGCCCTGCTCCAGGAGCACCTCCGAGAGCACAGCGGCCCTGGGCTGCCTGGTCAAGGACTACAAAA"
    
    find_best_motif(seq1)
    #inidicate_startpoint()

# perhaps make the output as texshade?




#if (__name__ == "__main__"):
#    status = main(sys.argv)
#    sys.exit(status)


