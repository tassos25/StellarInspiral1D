#!/bin/bash

#step1
cd step1_HMS-HMS
./clean
./mk
./rn | tee screen.txt
./cpy
cd ..

#step2
cd step2_RelaxModel
./clean
./mk
./rn | tee screen.txt
./cpy
cd ..

#step3
cd step3_inspiral
./clean
./mk
./rn | tee screen.txt
./cpy
cd ..
