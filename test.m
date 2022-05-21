% Just run this file as an example of WSCE

%% Load sample data set - The CNAE9 data set  
clear all;close all;clc;
load('SampleDataSet(CNAE9).mat');
% Data is unlabeled samples
% Class is ground truth (start from 1 to ...)
%% Run Weighted Spectral Cluster Ensemble Algorithm
% First method:
%Index = WSCE(Data, max(Class));
% Second method (recommmended)
block_size = 10;
num_neighbors = 10;
ShowDendrogram = 1;
Index = WSCE(Data, max(Class), block_size, num_neighbors, ShowDendrogram); 
%% Calculating the Classification Accuracy
ClassificationAccuracy = accuracy(Class, Index);
disp('The accuracy is: ');
disp(ClassificationAccuracy);