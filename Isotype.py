#!/usr/bin/python

"""Final Project for introscicomp2014: Antibody Isotype calculator.  Problem: our lab uses a very slow and inefficient program to determine antibody sequences.  Other labmates had developed means to analyze mutation and gene repertoire information but currently there is no way to properly determine the antibody isotype.  This program is designed to scan FASTA sequence files of immunologobulin genes and provide a best match, by percentage, to the antibody isotype that is likely being utilized. IgM, IgG1, IgG2, IgG3, IgG4, IgA1, IgA2, IgD, and IgE are all represented in our readout"""

__author__ = 'Noel Pauli (npauli@uchicago.edu)'
__version__ = '0.0.1'

#imports
"""import sys
import Bio
import Numpy"""

# this program will utilize Numpy and Biopython to do the intial sequence analysis

# constants

# functions

# These are the two sequences to match
"""seq2 = "CATCCCCGACCAGCCCCAAGGTCTTCCCGCTGAGCCTCGACA"
"""


seq1 = "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNATNCNGGTAGCTGGGGGGTCCTCAGGCTCTCCTGTGTAGTCTCTGGACTCACTTTCGACAATGCCTGGATGAGCTGGGTCCGCCAGGCTCCAGGGAAGGGGCTGGAGTGGGTTGGCCGTATTAAAAGCAAATCTGATGGTGGGACAACAGACTTCGCTGCACCCGTGAAAGGCAGATTCACCATCTCTAGAGATGACTCAAAAAATACAGTGTTTCTGCAAATGAACAGCCTGCAGACCGAGGACACAGCCGTGTATTACTGTGCCACAGCCCCCGGATTCCACTATTATGCTCCCTTTGACTACTGGGGCCCGGGAACCCTGGTCACCGTCTCCTCAGCCTCCACCAAGGGCCCATCGGTCTTCCCCCTGGCGCCCTGCTCCAGGAGCACCTCCGAGAGCACAGCGGCCCTGGGCTGCCTGGTCAAGGACTACAAAA"

IgG1 = "CCTCCACCAAGGGCCCATCGGTCTTCCCCCTGGCACCCTCCT"
IgG2 = "CCTCHACCAAGGGCCCATCGGTCTTCCCCCTGGCGCCCTGCT"
IgG3 = "CTTCCACCAAGGGCCCATCGGTCTTCCCCCTGGCGCCCTGCT"
IgG4 = "CTTCCACCAAGGGCCCATCCGTCTTCCCCCTGGCGCCCTGCT"
IgA1 = "CATCCCCGACCAGCCCCAAGGTCTTCCCGCTGAGCCTCTGCA"
IgM = "GGAGTGCATCCGCCCCAACCCTTTTCCCCCTCGTCTCCTGTG"
IgA2 = "CATCCCCGACCAGCCCCAAGGTCTTCCCGCTGAGCCTCGACA"
IgD = "CACCCACCAAGGCTCCGGATGTGTTCCCCATCATATCAGGGT"
IgE = "CCTCCACACAGAGCCCATCCGTCTTCCCCTTGACCCGCTGCT"
H = "HHHHHXHXHXXHX"


motifs = [IgG1, IgG2, IgG3, IgG4, IgA1, IgM, IgA2, IgD, IgE]
"""motifs = "H"
"""
# assign the longest sequence
# to s1, and the shortest to S2
# L1 is the length of the longest,
# l2 that of the shortest
"""for i in motifs:
    l2 = len(motifs)"""

l1 = len(seq1)
l2 = len(motifs[0])
#len simply counts the number of characters in the string

if l1 >= l2:
    s1 = seq1
    s2 = str(motifs)
# If the length of seq1 is greater than seq2, do basically nothing

else:
    s1 = seq2
    s2 = seq1
    l1, l2 = l2, l1 # swap the two lengths
# If the length of seq2 is greater than seq 1, swap the sequences

# function that computes a score
# by returning the number of matches
# starting from arbitrary startpoint
def calculate_score(s1, s2, l1, l2, startpoint):
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

# now try to find the best match (highest score)
my_best_align = None
my_best_score = -1
my_best_ratio = 0

for x in motifs:
    for i in range(l1):
        z = calculate_score(s1, x, l1, l2, i)
        print z
        if z [0] > my_best_score:
            my_best_align = "." * i + x
            my_best_score = z [0]
            my_best_ratio = z [1]
# z is a tuple in this case with two parameters. 
        elif z [0] == my_best_score: 
    	    if z [1] > my_best_ratio:
    	    	    my_best_align = "." * i + x
                    my_best_score = z [0]
                    my_best_ratio = z [1]   

print my_best_align
print s1
print "Best score:", my_best_score
print "Best Ratio:", my_best_ratio





"""    
if (__name__ == "__main__"):
    status = main(sys.argv)
    sys.exit(status)
"""
"""from Bio.Seq import Seq
    IgG1 = Seq("CCTCCACCAAGGGCCCATCGGTCTTCCCCCTGGCACCCTCCT")
    IgG2 = Seq("CCTCCACCAAGGGCCCATCGGTCTTCCCCCTGGCGCCCTGCT")
    IgG3 = Seq("CTTCCACCAAGGGCCCATCGGTCTTCCCCCTGGCGCCCTGCT")
    IgG4 = Seq("CTTCCACCAAGGGCCCATCCGTCTTCCCCCTGGCGCCCTGCT")
    IgA1 = Seq("CATCCCCGACCAGCCCCAAGGTCTTCCCGCTGAGCCTCTGCA")
    IgM = Seq("GGAGTGCATCCGCCCCAACCCTTTTCCCCCTCGTCTCCTGTG")
    IgA2 = Seq("CATCCCCGACCAGCCCCAAGGTCTTCCCGCTGAGCCTCGACA")
    IgD = Seq("CACCCACCAAGGCTCCGGATGTGTTCCCCATCATATCAGGGT")
    IgE = Seq("CCTCCACACAGAGCCCATCCGTCTTCCCCTTGACCCGCTGCT")"""
