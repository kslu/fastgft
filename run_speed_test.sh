#!/bin/bash

INPUT10="input/data10_20000.txt"
INPUT12="input/data12_20000.txt"
INPUT15="input/data15_20000.txt"
INPUT16="input/data16_20000.txt"
INPUT25="input/data25_20000.txt"
INPUT64="input/data64_20000.txt"
INPUT80="input/data80_20000.txt"
INPUT100="input/data100_20000.txt"

./speed_star10 $INPUT10 output/rt_star10_20000.txt
./speed_star100 $INPUT100 output/rt_star100_20000.txt
./speed_cycle12 $INPUT12 output/rt_cycle12_20000.txt
./speed_cycle80 $INPUT80 output/rt_cycle80_20000.txt
./speed_bd4x4 $INPUT16 output/rt_bd4x4_20000.txt
./speed_bd8x8 $INPUT64 output/rt_bd8x8_20000.txt
./speed_z4x4 $INPUT16 output/rt_z4x4_20000.txt
./speed_z8x8 $INPUT64 output/rt_z8x8_20000.txt
./speed_sk15 $INPUT15 output/rt_sk15_20000.txt
./speed_sk25 $INPUT25 output/rt_sk25_20000.txt
