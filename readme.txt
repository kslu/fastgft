Fast GFT implementation

This repository is provided for generating the results in the paper:

Keng-Shih Lu and Antonio Ortega, "Fast Graph Fourier Transforms based on
Graph Symmetry and Bipartition."

Here is a brief summary of how to produce the main results.

1. The C code is used to obtain the empirical run time of the fast and
   matrix GFTs. To build the binary files, type
	     
     make

	 and run the programs using commands like

	   ./speed_bd8x8 input/data64_20000.txt output/runtime_bd8x8_20000.txt

   where the first argument above specifies the input data file, and the
	 second argument specifies the file to be output. Input data of 20000
	 samples used in the paper is given in the input/ folder.

2. To generate the figures, use the MATLAB code in the matlab/ folder.
   The main script is demo_rtvserr.m

	 Before running this script, one can obtain the relative error (RE) and
	 number of multiplications of each GFT implementation purely in MATLAB.
	 To get these results, run the script demo_get_re_nmult.m before running
	 demo_rtvserr.m

3. All the GFT matrices and Givens rotations in the C code

   src/gftmatrices.h
	 src/givens_ptj.h
	 src/givens_tj.h

	 were generated using the MATLAB script demo_generate_c_code.m
	 (The resulting C code will be stored in matlab/c_code/ and should be
	 the same as those written the above .h files)

The programs and MATLAB scripts have been tested on macOS High Sierra 10.13.3
