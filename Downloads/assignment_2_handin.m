%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solution to Assignment, Project or Exam
% Populate Q structs with answers (Assig / Exam)
% Populate S structs with answers (Project)
% 
% Once you are done ...
% Run the solution script (this script).
% Run makeMat341.p to generate a MAT file.
%
% Always include a comment above a Matlab function.
% It is automatically turned into a help screen.
%
% Try it. To see this comment:
% In Matlab, enter: help solveAssig
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% A blank line terminates the help screen
% This comment does not show up when you type 'help'

% Clean up workspace
clear all; clc;

% SN variable must contain Student Number
% This must be right or solution will not be graded
SN    = 14425508;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% START OF FILE

s = tf('s')

% Q1
Q1.Tr = 0.136
Q1.Tp = 0.381
Q1.Ts = 0.807
Q1.Pos = 42.5871

% Q2
Q2.Z = 0.2622
Q2.Wn = 13990.3
Q2.G = (5.480*10^9/(s^2 + 7337*s+1.957*10^8))

% Q3
Q3.Wn = 8544.6
Q3.G = 2.044*10^9 / (s^2 + 4481*s + 7.301*10^7)

% Q4
Q4.Z = 0.3720
Q4.Wn = 8.8832*10^3
Q4.G = (2.21*10^9/(s^2 + 6609*s + 7.891*10^7))

% Q5
Q5.Z = 0.5278
Q5.Wn = 9.708*10^3
Q5.G = 2.639*10^9/(s^2 + 1.025*10^4*s + 9.425*10^7)

% Q6
Q6.Z = 1
Q6.Wn = 4956.6
Q6.G = (6.879*10^8/(s^2 + 9913*s + 2.457*10^7))

% Q7
Q7.Tr = 0.1795
Q7.Te = 0.0435
Q7.G = (3.149*10^9/(s^2 + 5559*s + 1.124*10^8))

%%% END OF FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
