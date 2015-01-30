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
instead of looping over n<sup>2</sup> or n*xm values
we loop only over



Results
================

Strikingly the K Banded algorithm does not show show improvement over 
needleman wunsch algorithm for the test cases we have in `tests/data`

The run time is depdicted by

![plot](./comparison.png)
