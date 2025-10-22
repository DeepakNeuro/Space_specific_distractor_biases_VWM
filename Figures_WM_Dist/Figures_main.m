%{
Main figure script for other figures and loading behavior data . 
%}

% clc; 
% clear all; 
%% path settings
addpath(genpath('/media/hdd/Sanchit/Exogenous_Project/Data Analysis Codes/EEG Analysis'));
addpath '/media/hdd/Sanchit/Exogenous_Project/Data Analysis Codes';
addpath(genpath("/media/hdd/Sanchit/Exogenous_Project/Toolboxes/MyFunctions"));

addpath(genpath('./MyFunctions'));
%% Load behavior data
beh = load('./Data/Behaviour/DataAnalysis_Bias_March_10_2023_11_23_08.mat');

nSubjects = length(beh.Subjects);
nEEG = 23; 
nSim = 10; 
[idxL, idxR] = beh_findAttendSide(beh.subjectNumber, beh.trialData); 
idxnoDist = beh.trialData.distractorSide==3;

%% Default plot settings
plotDefaults; 
tl = tiledlayout(8,6,'TileSpacing','compact', 'Padding', 'compact');


%% color maps for various figures
color_map = lines(5);
color_map(6,:) = [0 0 0]; 

color_map_att = brewermap(5, 'PiYG');
color_map_att = color_map_att([1, 5],:); 

color_map_split = color_map(3:4,:); 

% %% Figure Scripts
% F1; % Fig. 1. Behavioral plots; 
% F2; % Fig. 2. Decoding and Neural Bias Plots
% F3; % Fig. 3. Distractor decoding and ERPs
% F4; % Fig. 4. Correlation Plots
% F5and6A_F; % Fig. 5 and 6A-F simulations results
% F6G_I; % Fig. 6G-I simulations results


