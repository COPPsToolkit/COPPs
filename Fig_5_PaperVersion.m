% Script to create Figure 5 in 
% **de Groot, O., F. Mazelis, R. Motto, A. Ristiniemi**
% "A Toolkit for Computing Constrained Optimal Policy Projections (COPPs)"
%% Preamble
clear
addpath('Toolkit');

%% Run individual script
run('IndividualProjections\Run_COPPs_Fig_5_a')
projections_all.simple_rules = projections.simple_rules; 
projections_all.COPPs_5a = projections.COPPs; 
save('FigureInfo','projections_all')
%% Run individual script
run('IndividualProjections\Run_COPPs_Fig_5_b')
load('FigureInfo','projections_all')
projections_all.COPPs_5b = projections.COPPs; 
save('FigureInfo','projections_all')
%% Run individual script
run('IndividualProjections\Run_COPPs_Fig_5_c')
load('FigureInfo','projections_all')
projections_all.COPPs_5c = projections.COPPs; 
save('FigureInfo','projections_all')
%% Load relevant info
clear
load('FigureInfo')
%% Set preferences
params.plotting.Data = {projections_all.simple_rules.data   , '-'  , 1 ,[0,0,0], 'Baseline' ; ...
        projections_all.COPPs_5a.data          , '--' , 1 ,[0,0,1], 'H=1'; ...
        projections_all.COPPs_5b.data          , '-.' , 1 ,[0,1,1], 'H=2'; ...
        projections_all.COPPs_5c.data          , ':'  , 1 ,[0,1,0], 'H=8'; ...
        };

params.plotting.VarsToPlot = {...
    'pinf_ann',[-3,1]  ,'Inflation'  ;...
    'og'      ,[-1,3]  ,'Output gap' ;...
    'r_ann'   ,[-10,10],'Policy rate';...
    };

params.plotting.P_past   = 4;
params.plotting.P_future = 13;

%% Plot chart
PlotCOPPs(projections_all, params);