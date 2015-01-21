## Algorithms for Computational Biology

## Build Status [![Build Status](https://travis-ci.org/saketkc/comp-bio.svg)](https://travis-ci.org/saketkc/comp-bio)

## Implemented Algorithms:

- [Global Alignment](docs/GlobalAlignment.md)


## How to run?
```
git clone git@gitub.com:saketkc/comp-bio.git
cd comp-bio
make
make test
make install

./bin/global_alignment tests/data/align-1000bp-deletions.fasta
```

## Running tests
We use [Catch](https://github.com/philsquared/Catch) for running unit test cases.
```
make test
cd tests
./global_alignment_test 

===============================================================================
All tests passed (1 assertion in 1 test case)


```

## Memory Leaks? NO
```
$ valgrind ./bin/global_alignment tests/data/test.fasta 
==32079== Memcheck, a memory error detector
==32079== Copyright (C) 2002-2013, and GNU GPL’d, by Julian Seward et al.
==32079== Using Valgrind-3.10.0.SVN and LibVEX; rerun with -h for copyright info
==32079== Command: ./bin/global_alignment tests/data/test.fasta
==32079== 

--------------------------------------------

Sequence 1: GAATTCAGTTA
Sequence 2: GGATCGA

--------------------------------------------

Sequence2 Length: 11
Sequence2 Length: 7

----------------------Optimal Alignment Start--------------------------

GAATTCAGTTA

GGA-TC-G--A

----------------------Optimal Alignment End--------------------------
Score: 3
==32079== 
==32079== HEAP SUMMARY:
==32079==     in use at exit: 0 bytes in 0 blocks
==32079==   total heap usage: 54 allocs, 54 frees, 10,772 bytes allocated
==32079== 
==32079== All heap blocks were freed -- no leaks are possible
==32079== 
==32079== For counts of detected and suppressed errors, rerun with: -v
==32079== ERROR SUMMARY: 0 errors from 0 contexts (suppressed: 0 from 0)

```
## Third Party Libraries

|Library   	| What for?   	| License   	|
|---	|---	|---	|
|[SimpleIni](https://github.com/brofield/simpleini)   	|  Parsing `.ini` files 	| MIT   	|
|[Catch](https://github.com/philsquared/Catch)   	|  Unit Tests 	| Boost 1.0  	|





## License

[GPL 3.0](LICENSE)
