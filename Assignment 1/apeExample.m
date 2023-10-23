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
% This line does not show up in the help screen

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize

% Clean up workspace
clear all; clc;

% SN variable must contain Student Number
% This must be right or solution will not be graded
SN    = 90300815;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculations

% Q1: Scalars
% Do not round scalars
% Keep answers as accurate as possible
Q1.K1 = 123;            % integer
Q1.K2 = 1e2;            % integer with exponent
Q1.K3 = 1.23;           % floating point
Q1.K4 = 1.23e-6;        % floating point with exponent
Q1.K5 = 1+j*2;          % complex

% Q2: Vectors
% May be specified as row or column vectors
% Length must be right
Q2.Va = [1 2 3];                  % 1x3 row vec
Q2.Vb = [1 2 3]';                 % Use transpose operator for col vec
Q2.Vc = [1                        % 3x1 col vec
         2
         3];

% Q3: Matrices
% Both dimensions must be right
Q3.M1 = [1 2                      % 2x3 matrix
         3 4
         5 6];
Q3.M2 = [1 2;3 4];                % Use ';' to create new row
Q3.M3 = eye(4)*1+j;               % Matrix functions & complex ok
Q3.M4 = eye(4)*tf('s');           % Laplace operator is ok

% Q4: LTI Objects
% Created using tf(), zpk() or Laplace Operator tf('s')
% NEVER use chgTimeUnit() - leave in units of (s)
Q4.C  = tf(1, [2 3]);             % Use tf()
Q4.G  = zpk(4, [5 6], 7);         % Use zpk()
Q4.H  = 2*tf('s')/(3*tf('s')+4);  % Explicit using s operator

% Q5: Text Message
% Enclose in single quotes.
Q5.X1 = 'Advice is either useful or pleasant, but never both.';

% Use sprintf() to include line breaks and variables.
% See the help screen on 'Formatting Text' to include special characters
% like the single-quote symbol.
pos = 100;
neg = -pos;
Q5.X2 = sprintf('\nIf you think %d and %d are basically the same number,\nlet''s apply that logic to your bank account and see if you\nchange your mind.\n', pos, neg);

% Use the display() function to display formatted text.
display(Q5.X1);
display(Q5.X2);

% End of file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
