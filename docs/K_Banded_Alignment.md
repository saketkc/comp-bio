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
algorithm is similar
