#!/bin/bash

#step1
cd step1_HMS-CO
./clean
./mk
./rn
./cpy
cd ..

#step2
cd step2_RelaxModel
./clean
./mk
./rn
./cpy
cd ..
