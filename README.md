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
## Third Party Libraries

|Library   	| What for?   	| License   	|
|---	|---	|---	|
|[SimpleIni](https://github.com/brofield/simpleini)   	|  Parsing `.ini` files 	| MIT   	|
|[Catch](https://github.com/philsquared/Catch)   	|  Unit Tests 	| Boost 1.0  	|





## License

[GPL 3.0](LICENSE)
