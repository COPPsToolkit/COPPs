% Script to create Figure C3 in 
% **de Groot, O., F. Mazelis, R. Motto, A. Ristiniemi**
% "A Toolkit for Computing Constrained Optimal Policy Projections (COPPs)"
%% Preamble
clear
addpath('Toolkit');

%% Run individual script
run('IndividualProjections\Run_COPPs_Fig_C3_a')
projections_all.simple_rules = projections.simple_rules; 
projections_all.COPPs_C3a = projections.COPPs; 
save('FigureInfo','projections_all')
%% Run individual script
run('IndividualProjections\Run_COPPs_Fig_C3_b')
load('FigureInfo','projections_all')
projections_all.COPPs_C3b = projections.COPPs; 
save('FigureInfo','projections_all')
%% Load relevant info
clear
load('FigureInfo')
%% Set preferences
params.plotting.Data = {projections_all.simple_rules.data   , '-'  , 1 ,[0,0,0], 'Baseline' ; ...
        projections_all.COPPs_C3a.data          , '--'  , 1 ,[0,0,1],'Commitment: \alpha=0.6'; ...
        projections_all.COPPs_C3b.data          , '-.'  , 1 ,[0,1,1],'Small deviations: \alpha=0.6'; ...
        };

params.plotting.VarsToPlot = {...
    'pi_ann',[-40,5],'Inflation'  ;...
    'og'    ,[-30,5],'Output gap' ;...
    'r_ann' ,[-1,6] ,'Policy rate';...
    };


params.plotting.P_past   = 4;
params.plotting.P_future = 28;

%% Plot chart
PlotCOPPs(projections_all, params);