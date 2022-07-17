  _______   _______   _______   _______   
 |  ___  | |  ___  | |  ___  | |  ___  |  ______ 
 | |   |_| | |   | | | |   | | | |   | | |  ____|
 | |       | |   | | | |___| | | |___| | | |____ 
 | |    _  | |   | | |  _____| |  _____| |____  |
 | |___| | | |___| | | |       | |        ____| |
 |_______| |_______| |_|       |_|       |______|

A Toolkit for Computing Constrained Optimal Policy Projections
Version 1.1

de Groot, Mazelis, Motto, Ristiniemi
ECB Working Paper 2021 No. 2555
https://ideas.repec.org/p/ecb/ecbwps/20212555.html

Copyright (C) de Groot, Mazelis, Motto, Ristiniemi

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation,  version 3 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.

--------------------------------------------------------------

For further information, visit https://github.com/COPPsToolkit/COPPs

--------------------------------------------------------------

Content:


A. Files

1: 	License agreement, see LICENSE.txt

2: 	A manual that describes how to produce COPPs with new models and additional options, see COPPs_Manual_v1p1.pdf. 

3:	Run_COPPS.m can be copied to be used as a template. It provides a simple example that runs optimal policy projections for the Smets and Wouters model without constraints. 

4: 	Replication files for all figures in the working paper. 
	Open the script Fig_XY_PaperVersion.m related to chart XY and execute. The Matlab paths to the optimization toolbox and Dynare need to be set manually. 

5: 	Replication file for Table D.1 in the working paper on computational times. 
	Open the script Tab_1.m and execute.

6: 	Change log, see CHANGELOG.txt


B. Folders

1:	\Toolkit includes the core files required to execute any Constrained Optimal Policy Projection.

2:	\Models contains several folders with dynare mod-files and related datafiles, which are drawn on in the computation of COPPs. 

3:	\IndividualProjections includes all individual COPPs scripts that are drawn on for the replication of figures in the working paper. 

4:	\MiscReplication provides the mod-files that compute optimal policy via the Lagrangian approach. 

