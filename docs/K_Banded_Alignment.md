K Band Global Alignment
=============================


K Band global alignment uses a clever approach
of calculating the alignment scores for elements
close to the diagonal. The underlying hypothesis
being that an alignment of two sequences
with at most _k_ differences(mismatches+indels)
will always stay close to the diagonal.

Approach
=================

The underlying implementation for K banded
algorithm is similar to Needleman Wunsch algorithm for global
alignment as detailed in [GlobalAlignment.md](./GlobalAlignment.md) but
instead of looping over n<sup>2</sup> or n\*m values
we loop only over


Implementation
=====================
```
D[i,0] = i\*INDEL

D[0,j] = j\*INDEL

score  = infinity(some large value)

minDistance = 0

k = 1

while(score > minDistance){

    performKBandedAlignment()

    k*=2

    minDistance = (seq2.length()-k-1)*dMATCH + (2*(k+1)+abs(seq1.length()-seq2.length()))*dINDEL;

    score = getOptimalScore(SM);


}

getoptimalscore() returns D[seq1_length, seq2_length] and seq1_length >= seq2_length is ensured.

```

Expected running time
==============================

`kn + (1\2)kn + (1\4)kn + ... <=2kn`

So ~O(2n<sup>2</sup>)

Results
================
Not upto the mark.

To replicated the results: Run `bash benchmark.sh`

Strikingly the K Banded algorithm does not show show improvement over 
needleman wunsch algorithm for the test cases we have in `tests/data`

The run time is depicted by

![plot](./comparison.png)

