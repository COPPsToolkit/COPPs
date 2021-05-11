% Script to create Figure C1 in 
% **de Groot, O., F. Mazelis, R. Motto, A. Ristiniemi**
% "A Toolkit for Computing Constrained Optimal Policy Projections (COPPs)"
%% Preamble
clear
addpath('Toolkit');

%% Run individual script
run('IndividualProjections\Run_COPPs_Fig_C1_a')
projections_all.simple_rules = projections.simple_rules; 
projections_all.COPPs_C1a = projections.COPPs; 
save('FigureInfo','projections_all')
%% Run individual script
run('IndividualProjections\Run_COPPs_Fig_C1_b')
load('FigureInfo','projections_all')
projections_all.COPPs_C1b = projections.COPPs; 
save('FigureInfo','projections_all')
%% Load relevant info
clear
load('FigureInfo')
%% Set preferences
params.plotting.Data = {projections_all.simple_rules.data   , '-'  , 1 ,[0,0,0], 'Baseline' ; ...
        projections_all.COPPs_C1a.data          , '--'  , 1 ,[0,0,1],'Commitment'; ...
        projections_all.COPPs_C1b.data          , '-.' , 1 ,[0,1,0], 'Discretion'; ...
        };

params.plotting.VarsToPlot = {...
    'pinf_ann',[-2,4] ,'Inflation'  ;...
    'og'      ,[-4,2] ,'Output gap' ;...
    'r_ann'   ,[-2,12],'Policy rate';...
    };

params.plotting.P_past   = 4;
params.plotting.P_future = 9;

%% Plot chart
PlotCOPPs(projections_all, params);