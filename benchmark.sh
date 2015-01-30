#!/bin/bash

for file in `ls -1 tests/data/`; do
   cmd="./bin/global_alignment_distancewise tests/data/${file} >> global.csv"
   eval $cmd
   eval "./bin/banded_global_alignment tests/data/$file >> banded.csv"
   done

