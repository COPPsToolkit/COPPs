% Script to calculate timings
% **de Groot, O., F. Mazelis, R. Motto, A. Ristiniemi**
% "A Toolkit for Computing Constrained Optimal Policy Projections (COPPs)"
%% Preamble
clear

tic;Run_COPPs_Fig_2_a;toc;close all;clc
proj.fig2a.params   = params;
proj.fig2a.toc      = toc;
save('TimingInfo','proj')

tic;Run_COPPs_Fig_2_b;toc;close all;clc
load('TimingInfo','proj')
proj.fig2b.params   = params;
proj.fig2b.toc      = toc;
save('TimingInfo','proj')

tic;Run_COPPs_Fig_4_a;toc;close all;clc
load('TimingInfo','proj')
proj.fig4a.params   = params;
proj.fig4a.toc      = toc;
save('TimingInfo','proj')

tic;Run_COPPs_Fig_4_b;toc;close all;clc
load('TimingInfo','proj')
proj.fig4b.params   = params;
proj.fig4b.toc      = toc;
save('TimingInfo','proj')

tic;Run_COPPs_Fig_8_b;toc;close all;clc
load('TimingInfo','proj')
proj.fig8b.params   = params;
proj.fig8b.toc      = toc;
save('TimingInfo','proj')

tic;Run_COPPs_Fig_8_c;toc;close all;clc
load('TimingInfo','proj')
proj.fig8c.params   = params;
proj.fig8c.toc      = toc;
save('TimingInfo','proj')

tic;Run_COPPs_Fig_9_a;toc;close all;clc
load('TimingInfo','proj')
proj.fig9a.params   = params;
proj.fig9a.toc      = toc;
save('TimingInfo','proj')

tic;Run_COPPs_Fig_10_a;toc;close all;clc
load('TimingInfo','proj')
proj.fig10a.params   = params;
proj.fig10a.toc      = toc;
save('TimingInfo','proj')

tic;Run_COPPs_Fig_11_a;toc;close all;clc
load('TimingInfo','proj')
proj.fig11a.params   = params;
proj.fig11a.toc      = toc;
save('TimingInfo','proj')

tic;Run_COPPs_Fig_11_c;toc;close all;clc
load('TimingInfo','proj')
proj.fig11c.params   = params;
proj.fig11c.toc      = toc;
save('TimingInfo','proj')

tic;Run_COPPs_Fig_13_a;toc;close all;clc
load('TimingInfo','proj')
proj.fig13a.params   = params;
proj.fig13a.toc      = toc;
save('TimingInfo','proj')

tic;Run_COPPs_Fig_13_b;toc;close all;clc
load('TimingInfo','proj')
proj.fig13b.params   = params;
proj.fig13b.toc      = toc;
save('TimingInfo','proj')

%%
tic;Run_COPPs_Fig_4_a_tempH20;toc;close all;clc
load('TimingInfo','proj')
proj.fig4aH20.params   = params;
proj.fig4aH20.toc      = toc;
save('TimingInfo','proj')

tic;Run_COPPs_Fig_4_a_tempH40;toc;close all;clc
load('TimingInfo','proj')
proj.fig4aH40.params   = params;
proj.fig4aH40.toc      = toc;
save('TimingInfo','proj')

tic;Run_COPPs_Fig_4_b_tempH20;toc;close all;clc
load('TimingInfo','proj')
proj.fig4bH20.params   = params;
proj.fig4bH20.toc      = toc;
save('TimingInfo','proj')

tic;Run_COPPs_Fig_4_b_tempH40;toc;close all;clc
load('TimingInfo','proj')
proj.fig4bH40.params   = params;
proj.fig4bH40.toc      = toc;
save('TimingInfo','proj')

%%
fprintf(' 2 & NK3 & 1 & COM & UNC & %2.0f & %2.2f\\\\\n',[proj.fig2a.params.H_policy,proj.fig2a.toc]);
fprintf(' 2 & NK3 & 1 & DIS & UNC & %2.0f & %2.2f\\\\\n',[proj.fig2b.params.H_policy,proj.fig2b.toc]);
fprintf(' 4 & NK3 & 1 & COM & CON & %2.0f & %2.2f\\\\\n',[proj.fig4aH20.params.H_policy,proj.fig4aH20.toc]);
fprintf(' 4 & NK3 & 1 & COM & CON & %2.0f & %2.2f\\\\\n',[proj.fig4aH40.params.H_policy,proj.fig4aH40.toc]);
fprintf(' 4 & NK3 & 1 & COM & CON & %2.0f & %2.2f\\\\\n',[proj.fig4a.params.H_policy,proj.fig4a.toc]);
fprintf(' 4 & NK3 & 1 & DIS & CON & %2.0f & %2.2f\\\\\n',[proj.fig4bH20.params.H_policy,proj.fig4bH20.toc]);
fprintf(' 4 & NK3 & 1 & DIS & CON & %2.0f & %2.2f\\\\\n',[proj.fig4bH40.params.H_policy,proj.fig4bH40.toc]);
fprintf(' 4 & NK3 & 1 & DIS & CON & %2.0f & %2.2f\\\\\n',[proj.fig4b.params.H_policy,proj.fig4b.toc]);
fprintf(' 8 & NK3 & 1 & LC (D=6) & CON & %2.0f & %2.2f\\\\\n',[proj.fig8b.params.H_policy,proj.fig8b.toc]);
fprintf(' 8 & NK3 & 1 & LC (D=9) & CON & %2.0f & %2.2f\\\\\n',[proj.fig8c.params.H_policy,proj.fig8c.toc]);
fprintf(' 9 & SW07 & 1 & COM & UNC & %2.0f & %2.2f\\\\\n',[proj.fig9a.params.H_policy,proj.fig9a.toc]);
fprintf(' 10 & SW07 & 1 & COM & CON & %2.0f & %2.2f\\\\\n',[proj.fig10a.params.H_policy,proj.fig10a.toc]);
fprintf(' 11 & SWu20 & 1 & COM & CON & %2.0f & %2.2f\\\\\n',[proj.fig11a.params.H_policy,proj.fig11a.toc]);
fprintf(' 11 & SWu20 & 2 & COM & CON & %2.0f & %2.2f\\\\\n',[proj.fig11c.params.H_policy,proj.fig11c.toc]);
fprintf(' 13 & SWu20 & 1 & DIS & CON & %2.0f & %2.2f\\\\\n',[proj.fig13a.params.H_policy,proj.fig13a.toc]);
fprintf(' 13 & SWu20 & 2 & DIS & CON & %2.0f & %2.2f\\\\\n',[proj.fig13b.params.H_policy,proj.fig13b.toc]);

