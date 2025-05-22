# pdtrajectories
This repository contains code used to create fitted trajectory plots based on 
linear mixed effects models for a project on Parkinson's disease.

Input options include the total data set (multiple measurements per person),
a data set with just one measurement per person with baseline covariates,
a binary "split" variable to display two different trajectories on the graphs,
a customized figure title, and the outcome variable.

Much of this code is written specifically for Parkinson's cognitive models for
an RA project in 2025. Model predictors can be adjusted for other projects as
desired, but currently are not function inputs 
(the actual function code should be modified to do this.)

This code was tested using R version 4.4.2.
