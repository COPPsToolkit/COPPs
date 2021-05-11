% Script to create Figure C2 in 
% **de Groot, O., F. Mazelis, R. Motto, A. Ristiniemi**
% "A Toolkit for Computing Constrained Optimal Policy Projections (COPPs)"
%% Preamble
clear
addpath('Toolkit');

%% Run individual script
run('IndividualProjections\Run_COPPs_Fig_C2_a')
projections_all.simple_rules = projections.simple_rules; 
projections_all.COPPs_C2a = projections.COPPs; 
save('FigureInfo','projections_all')
%% Load relevant info
clear
load('FigureInfo')
%% Set preferences
params.plotting.Data = {projections_all.simple_rules.data   , '-'  , 1 ,[0,0,0], 'Baseline' ; ...
        projections_all.COPPs_C2a.data          , '--'  , 1 ,[0,0,1],'Optimal policy projection'; ...
        };

params.plotting.VarsToPlot = {...
    'pi_p_ann',[-.8,.2] ,'Inflation'  ; ...
    'y_gap'   ,[-.6,.05],'Output gap' ; ...
    'r'       ,[-.3,.1] ,'Policy rate'; ...
    };

params.plotting.P_past   = 4;
params.plotting.P_future = 28;

%% Plot chart
PlotCOPPs(projections_all, params);