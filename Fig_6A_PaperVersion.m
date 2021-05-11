% Script to create Figure 6A in 
% **de Groot, O., F. Mazelis, R. Motto, A. Ristiniemi**
% "A Toolkit for Computing Constrained Optimal Policy Projections (COPPs)"
%% Preamble
clear
addpath('Toolkit');

%% Run individual script
run('IndividualProjections\Run_COPPs_Fig_6_a')
projections_all.simple_rules = projections.simple_rules; 
projections_all.COPPs_6a = projections.COPPs; 
save('FigureInfo','projections_all')
%% Run individual script
run('IndividualProjections\Run_COPPs_Fig_6_b')
load('FigureInfo','projections_all')
projections_all.COPPs_6b = projections.COPPs; 
save('FigureInfo','projections_all')
%% Run individual script
run('IndividualProjections\Run_COPPs_Fig_6_c')
load('FigureInfo','projections_all')
projections_all.COPPs_6c = projections.COPPs; 
save('FigureInfo','projections_all')
%% Load relevant info
clear
load('FigureInfo')
%% Set preferences
params.plotting.Data = {projections_all.simple_rules.data   , '-'  , 1 ,[0,0,0], 'Baseline' ; ...
        projections_all.COPPs_6a.data          , '--' , 1 ,[0,0,1], '\alpha=1'; ...
        projections_all.COPPs_6b.data          , '-.' , 1 ,[0,1,1], '\alpha=0.6'; ...
        projections_all.COPPs_6c.data          , ':'  , 1 ,[0,1,0], '\alpha=0'; ...
        };

params.plotting.VarsToPlot = {...
    'pi_ann',[-40,5],'Inflation'  ;...
    'og'    ,[-30,5],'Output gap' ;...
    'r_ann' ,[-1,6] ,'Policy rate';...
    };

params.plotting.P_past   = 4;
params.plotting.P_future = 25;

%% Plot chart
PlotCOPPs(projections_all, params);