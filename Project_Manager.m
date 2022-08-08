%% Project Note


%%

clear; 
clc; 
close all;


% Include dependencies
addpath('./calibration data'); 
addpath(genpath('./packages'));
addpath('./common'); 
addpath('./make odes'); 
addpath(genpath('./shared results')); 
addpath('./user files'); 

rootwd = pwd;


% Select a simulation tas from the list

listJob = {
    'model calibration';        % (1)
    'drug response simulation';  % (2)
    'calculate synergy score';   % (3)
    };

[ jobID ] = readInput( listJob );
jobcode = listJob{jobID}; % Selected

Task_Manager
