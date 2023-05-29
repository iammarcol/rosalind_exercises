# Given two strings s s  and t t  of equal length, return the Hamming distance between s s  and t t , 
# denoted dH(s,t) d H ( s , t ) , is the number of corresponding symbols that differ in s s  and t t 
# input:GAGCCTACTAACGGGAT
#       CATCGTAATGACGGCCT
def hamming_distance(fastafile):
    count=0
    sequence=[]
    with open(fastafile,'r') as file:
        for line in file:
            sequence.append(line.strip())
        for i in range(len(sequence[0])):
            aa1=sequence[0][i]
            for j in range(1,len(sequence)):
                aa2=sequence[j][i]
                if aa1!=aa2:
                    count+=1
        print(count)

# hamming_distance("inputex1.fasta")

###################################################################################################

# For DNA strings s1 s 1  and s2 s 2  having the same length, their transition/transversion ratio R(s1,s2) R ( s 1 , s 2 )  
# is the ratio of the total number of transitions to the total number of transversions, where symbol substitutions are inferred 
# from mismatched corresponding symbols as when calculating Hamming distance (see “Counting Point Mutations”).
# TRANSITIONS: A to G or G to A, C to T or T to C
# TRANSVERSIONS:  A to C, G to T, C to G, A to T and vice versa all

# count transitions, count transversions and then return the ratio

def mutation_ratio(fastafile):
    transitions = 0
    transversions = 0
    tt_ratio = 0
    sequence = ""
    seq_list = []
    with open(fastafile, 'r') as file:
        for line in file:
            if line.startswith(">"):
                if sequence:
                    seq_list.append(sequence.strip())
                identif = line.strip()
                sequence = ""
            else:
                sequence += line.strip()        
        if sequence:
            seq_list.append(sequence.strip())
    for i in range(len(seq_list[0])):
        amino_acids = set(seq[i] for seq in seq_list)
        if len(amino_acids) == 2:
            if {"A", "G"} <= amino_acids or {"C", "T"} <= amino_acids:
                transitions += 1
        if len(amino_acids) == 2 or len(amino_acids) == 3:
            if {"A", "C"} <= amino_acids or {"G", "T"} <= amino_acids or {"C", "G"} <= amino_acids or {"A", "T"} <= amino_acids:
                transversions += 1
    print(transitions)
    print(transversions)
    print(transitions/transversions)

# mutation_ratio("inputex2.fasta")

###################################################################################################

# For two strings s1 s 1  and s2 s 2  of equal length, the p-distance between them, denoted dp(s1,s2) d p ( s 1 , s 2 ) , 
# is the proportion of corresponding symbols that differ between s1 s 1  and s2 s 2 .  For a general distance function d d 
#  on n n  taxa s1,s2,…,sn s 1 , s 2 , … , s n  (taxa are often represented by genetic strings), we may encode the distances 
# between pairs of taxa via a distance matrix D D  in which Di,j=d(si,sj) D i , j = d ( s i , s j ) .  Given: A collection of 
# n n  (n≤10 n ≤ 10 ) DNA strings s1,…,sn s 1 , … , s n  of equal length (at most 1 kbp). Strings are given in FASTA format.  
# Return: The matrix D D  corresponding to the p-distance dp d p  on the given strings

# put all sequences in a seq_list, compare one by one seq iterating through the list

import numpy as np
def p_dist(fastafile):
    sequence = ""
    seq_list = []
    differences=0
    diff_list=[]
    with open(fastafile, 'r') as file:
        for line in file:
            if line.startswith(">"):
                if sequence:
                    seq_list.append(sequence.strip())
                identif = line.strip()
                sequence = ""
            else:
                sequence += line.strip()        
        if sequence:
            seq_list.append(sequence.strip())
        size = int(len(seq_list))
        for i in range(len(seq_list)):
            for j in range(len(seq_list)):
                aa1=seq_list[i]
                aa2=seq_list[j]
                common_letters = float((sum(1 for a, b in zip(aa1, aa2) if a != b))/len(aa1))
                diff_list.append(common_letters)
    matrix = np.reshape(diff_list, (size, size))
    print(matrix)

# p_dist("inputex3.fasta")

###################################################################################################


# Given two strings s s  and t t  (of possibly different lengths), the edit distance dE(s,t) d E ( s , t )  is 
# the minimum number of edit operations needed to transform s s  into t t , where an edit operation is defined as 
# the substitution, insertion, or deletion of a single symbol.  The latter two operations incorporate the case in 
# which a contiguous interval is inserted into or deleted from a string; such an interval is called a gap. For the 
# purposes of this problem, the insertion or deletion of a gap of length k k  still counts as k k  distinct edit operations.  
# Given: Two protein strings s s  and t t  in FASTA format (each of length at most 1000 aa)
# Return: The edit distance dE(s,t) d E ( s , t ) 

def edit_dist(fastafile):
    sequence = ""
    seq_list = []
    with open(fastafile, 'r') as file:
        for line in file:
            if line.startswith(">"):
                if sequence:
                    seq_list.append(sequence.strip())
                identif = line.strip()
                sequence = ""
            else:
                sequence += line.strip()        
        if sequence:
            seq_list.append(sequence.strip())
        s = seq_list[0]
        t = seq_list[1]
        rows = len(s)+1
        columns = len(t)+1
        dp = [[0] * columns for _ in range(rows)]
# fill in the matrix using dynamic programming
        for i in range(1, rows):
            dp[i][0] = i
        for j in range(1, columns):
            dp[0][j] = j
        for i in range(1, rows):
            for j in range(1,columns):
                if s[i-1]==t[j-1]:
                    dp[i][j]=dp[i-1][j-1]
                else:
                    substitution = dp[i - 1][j - 1] + 1
                    insertion = dp[i][j - 1] + 1
                    deletion = dp[i - 1][j] + 1
                    dp[i][j] = min(substitution, insertion, deletion)
        edit_distance=dp[len(s)][len(t)]
        return edit_distance

# edit_dist("inputex4.fasta") 

###################################################################################################

