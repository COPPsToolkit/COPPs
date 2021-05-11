%% ########################################################################
%% ################ Constrained Optimal Policy Projections ################
%% ########################################################################
%
% **de Groot, O., F. Mazelis, R. Motto, A. Ristiniemi**
% "A Toolkit for Computing Constrained Optimal Policy Projections (COPPs)"
% This  version: 1.0
%__________________________________________________________________________
%% ########################################################################
%% ##################             Preamble               ##################
%% ########################################################################
%% 
clear;
clc;

% Add toolkit to matlab path 
addpath('..\Toolkit');

%% ########################################################################
%% ##################        Required user input         ##################
%% ########################################################################
%% 

% Model details
params.ModelDirectory         = '..\Models\GaliNKStickW';
params.ModelFilename          = 'GaliNKStickW';
params.PlannerDiscountFactor  = 'betta';

% Target variables and policy variables
params.LossFunctionVariables       = {'pi_p';'pi_w';'y_gap'};
params.PolicyInstrumentsAndShocks  = {'r','er'};

% Timeseries
params.PastPeriods                 = 50;  

%% ########################################################################
%% ##################        Optional user input         ##################
%% ########################################################################
%%
% Forward guidance
params.alpha = 1;

params.T_loss   = 60;
params.H_policy = 60;

params.Constraint     = 'OFF';
params.RateMin        = 0;        % ELB in annual policy rate
params.PolicyRateStSt = 0;

siggma = 1;
varphi = 5;
phi_pi = 1.5;
phi_y  = 0.125;
theta_p= 3/4;
rho_nu = 0.5;
rho_z  = 0.5;
rho_a  = 0.9;
betta  = 0.99;
eta    = 3.77; %footnote 11, p. 115
alppha = 1/4;
epsilon_p = 9;
epsilon_w = 4.5;
theta_w   = 3/4;

lambda_p = 0.0215;
lambda_w = 0.0037;

W.y_gap = siggma+(varphi+alppha)/(1-alppha);
W.pi_p  = epsilon_p/lambda_p;
W.pi_w  = epsilon_w*(1-alppha)/lambda_w;


params.LossFunctionWeights = [W.pi_p,W.pi_w,W.y_gap];

params.plotting.P_future = 17; 

params.plotting.NoPlot     = 1;
params.plotting.VarsToPlot = {...
    'pi_p_ann',[-.8,.2] ,'Inflation'  ; ...
    'y_gap'   ,[-.6,.05],'Output gap' ; ...
    'r'       ,[-.3,.1] ,'Policy rate'; ...
    };

%% ########################################################################
%% #########        No user input required below this line        #########
%% ########################################################################
%% Set non-user defined values to default
params = SetRemainingParameters(params);

%% Solve for the policy projections
projections = SolveCOPPs(params);

%% Plot everything
params.plotting.P_past   = 4;
params.plotting.P_future = 28;
PlotCOPPs(projections, params);
