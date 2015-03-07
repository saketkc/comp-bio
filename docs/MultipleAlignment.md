Multiple Global Alignment
=============================
The approach uses an approxmiate algorithm to perform Multiple sequence alignment
by sequences are iteratively aligned along a tree. The tree is constructed
using hiererchial clustering.

Starting with a sequence profile for each sequence
where the rows consist of the nucleotides or aminoacids
and the columns correspond to the position in sequence. Each 
entry of this matrix represents the frequency of that base
at that position in the sequence.

Using profiles from each sequence a distance matrix is created
reflecting pairwise distance between each pair of sequence.

The sequences are merged starting from the closest sequences first and an alignment
is created for these merged sequences.



Approach
=================

Given P = p1, p2, .., pn and Q = q1, q2, .., qn where P and Q represent sequence profile 
the score using profile alignment is given by g(pi,qj) = \sigma s[k,l] * p_i[k] q_j[l] 
where k,l = set of nucleotides A,T,G,C or a ‘-’ 
and s[k,l] represents the score for aligning k, against l

The optimal global alignment between profiles using a distance metric is given by
D[i,j] = min(D[i-i,j]+g[pi,’-’], D[i,j-1]+g[‘-’,qj], D[i,j]+g[pi,qj])


Expected running time
==============================

For k sequences with largest sequence length n, the run time is O(kn^2)

