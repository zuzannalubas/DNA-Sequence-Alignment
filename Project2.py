#!/usr/bin/env python
# coding: utf-8

# In[7]:


import numpy as np

def needleman_wunsch(seq1, seq2, match=1, mismatch=-1, gap=-1):
    """Implementation of the Needleman-Wunsch algorithm for global DNA sequence alignment."""
    n, m = len(seq1), len(seq2)
    score_matrix = np.zeros((n+1, m+1))
    
    # Initialize the matrix
    for i in range(n+1):
        score_matrix[i][0] = i * gap
    for j in range(m+1):
        score_matrix[0][j] = j * gap
    
    # Fill the matrix
    for i in range(1, n+1):
        for j in range(1, m+1):
            match_score = score_matrix[i-1][j-1] + (match if seq1[i-1] == seq2[j-1] else mismatch)
            delete = score_matrix[i-1][j] + gap
            insert = score_matrix[i][j-1] + gap
            score_matrix[i][j] = max(match_score, delete, insert)
    
    # Traceback to obtain the alignment
    align1, align2 = "", ""
    i, j = n, m
    while i > 0 and j > 0:
        if score_matrix[i][j] == score_matrix[i-1][j-1] + (match if seq1[i-1] == seq2[j-1] else mismatch):
            align1 += seq1[i-1]
            align2 += seq2[j-1]
            i -= 1
            j -= 1
        elif score_matrix[i][j] == score_matrix[i-1][j] + gap:
            align1 += seq1[i-1]
            align2 += "-"
            i -= 1
        else:
            align1 += "-"
            align2 += seq2[j-1]
            j -= 1
    
    return align1[::-1], align2[::-1]

# Example DNA sequences
seq1 = "GATTACA"
seq2 = "GCATGCU"

alignment1, alignment2 = needleman_wunsch(seq1, seq2)
print("Global Alignment:")
print(alignment1)
print(alignment2)


# In[ ]:




