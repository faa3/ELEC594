# ELEC594

This repository contains the source files to implement Wavelet delineation to segment ECG recordings.

## DWTDelineatoin.m
The source code is contained in DWTDelineation.m as it contains the proper functions that receive ECG recordings as an input and then compute the R-peaks, and the T wave. In the end, the functions return the segementation points of the heartbeat recordings. In the end, the heartbeats are processed as probability distributions in a 1D array.

## MSDVisualDWT.m
This file will receive a selected set of patients' heartbeats and their respective computed Euclidean and Wasserstein distance in the format of matrices array. The result will be a visualization of the Multidimensional scaling embedding of the distance matrices.

## Latex files
The repo also includes the pdf of the final report detailing the theory, algorithm, and results. Additionally, there is also the source files in latex format.
