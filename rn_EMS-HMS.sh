#!/bin/bash

#step1
cd step1_build_EMS
./clean
./mk
./rn | tee screen.txt
./cpy
cd ..

#step2.1
cd step2.1_RelaxModel_EMS_hydro
./clean
./mk
./rn | tee screen.txt
./cpy
cd ..

#step2.2
cd step2.2_RelaxModel_EMS_hydro
./clean
./mk
./rn | tee screen.txt
./cpy
cd ..

#step3
cd step3_Inspiral_EMS
./clean
./mk
./rn | tee screen.txt
./cpy
cd ..
