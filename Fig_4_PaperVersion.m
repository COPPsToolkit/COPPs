% Script to create Figure 4 in 
% **de Groot, O., F. Mazelis, R. Motto, A. Ristiniemi**
% "A Toolkit for Computing Constrained Optimal Policy Projections (COPPs)"
%% Preamble
clear
addpath('Toolkit');

%% Run individual script
run('IndividualProjections\Run_COPPs_Fig_4_a')
projections_all.simple_rules = projections.simple_rules; 
projections_all.COPPs_4a = projections.COPPs; 
save('FigureInfo','projections_all')
%% Run individual script
run('IndividualProjections\Run_COPPs_Fig_4_b')
load('FigureInfo','projections_all')
projections_all.COPPs_4b = projections.COPPs; 
save('FigureInfo','projections_all')
%% Load relevant info
clear
load('FigureInfo')
%% Set preferences
params.plotting.Data = {projections_all.simple_rules.data   , '-'  , 1 ,[0,0,0], 'Baseline' ; ...
        projections_all.COPPs_4a.data          , '--' , 1 ,[0,0,1], 'Commitment'; ...
        projections_all.COPPs_4b.data          , '-.' , 1 ,[0,1,0], 'Discretion'; ...
        };

params.plotting.VarsToPlot = {...
    'pi_ann',[-40,4],'Inflation'  ;...
    'og'    ,[-30,2],'Output gap' ;...
    'r_ann' ,[-1,5] ,'Policy rate';...
    };

params.plotting.P_past   = 4;
params.plotting.P_future = 16;

%% Plot chart
PlotCOPPs(projections_all, params);