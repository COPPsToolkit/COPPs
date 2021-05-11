% Script to create Figure 11 in 
% **de Groot, O., F. Mazelis, R. Motto, A. Ristiniemi**
% "A Toolkit for Computing Constrained Optimal Policy Projections (COPPs)"
%% Preamble
clear
addpath('Toolkit');

%% Run individual script
run('IndividualProjections\Run_COPPs_Fig_11_a')
projections_all.simple_rules = projections.simple_rules; 
projections_all.COPPs_11a = projections.COPPs; 
save('FigureInfo','projections_all')
%% Run individual script
run('IndividualProjections\Run_COPPs_Fig_11_b')
load('FigureInfo','projections_all')
projections_all.COPPs_11b = projections.COPPs; 
save('FigureInfo','projections_all')
%% Run individual script
run('IndividualProjections\Run_COPPs_Fig_11_c')
load('FigureInfo','projections_all')
projections_all.COPPs_11c = projections.COPPs; 
save('FigureInfo','projections_all')
%% Load relevant info
clear
load('FigureInfo')
%% Set preferences
params.plotting.Data = {projections_all.simple_rules.data   , '-'  , 1 ,[0,0,0], 'Baseline' ; ...
        projections_all.COPPs_11a.data          , ':' , 1 ,[0,0,1], 'SWu20: rates only'; ...
        projections_all.COPPs_11b.data          , '-.' , 1 ,[1,0,0], 'SWu20: qe only'; ...
        projections_all.COPPs_11c.data          , '--' , 1 ,[0,1,0], 'SWu20: rates+qe'; ...
        };

params.plotting.VarsToPlot = {...
    'obs_pinf_4q' , [ -1 ,  4 ]  , 'Inflation (annual, P.P.)'     ; ...
    'ygap'        , [ -9 ,  6 ]  , 'Output gap (P.P.)'            ; ...
    'obs_r_ann'   , [ -8 ,  7]   , 'Interest rate (annual, P.P.)' ; ...
    'qeobs'       , [ -5 ,  50 ] , 'Asset holdings (% of GDP) '   ; ...
    };

params.PastPeriods          = 219;  % 2009 March
params.T_full               = 78; % How far into the future to roll the forecast forward.
params.plotting.first_date  = datetime(1954,7,01); % Set the first date of the series
params.plotting.freq        = 'quarter'; % Frequency of data (month/quarter/year)

%% Plot chart
PlotCOPPs(projections_all, params);