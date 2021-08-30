This zip folder contains a number of generated data files, some real data, as 
well as the code used to produce or analyze the data herein. 

***R Code:
This code is used to produce any of the graphics used in the paper: 
Comparing Groups of Time Dependent Data Using Locally Weighted Scatterplot 
Smoothing Alpha-Adjusted Serial T-tests, by Tim Niiler, copyright 2019

Briefly:
buildDistributions.R builds distributions for false discovery rate testing.
buildFiguresFinal.R can be run to build any of the figures in the paper.
buildFiguresExtra.R can be run to build certain extra graphics that can 
    be considered for extra examples of the process.
    
Detailed commentary for all code can be found within the *.R files. 

***Python Code:
This code is used for RFT testing and is briefly described below:

testRFT.py - used to run RFT tests on data for direct comparison to LAAST
            in particular with the data used by Pataky et al.
nullTestRFT.py - used for FDR testing using RFT 
nullTestRFTskew.py - used for FDR testing on non-normal distributions using RFT

***Other Files
patakyT.txt - t-value data scanned from Pataky's figure comparing gastroc forces.
KinematicsEMGmomentarmForces-walking.csv - data from Bessier et al, and used
            by Pataky et al.  This is the basis for Figure 2 of the paper.

