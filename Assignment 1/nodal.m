%%%%%%%%%%%%%%%%%%%%
% Nodal Analysis (Matrix Method) Tutorial
%
% This is a stand-alone script for solving the RLC circuit problem
% that is used to introduce the Matrix Method for solving 
% nodal analysis problems in ELEC 341.
%
% Calling Syntax:
% nodal
%
% Note: All variables cleared when this is run
%
% Author: Leo Stocco
%%%%%%%%%%%%%%%%%%%%

clear all; clc;

% Parameters
R1 = 200;       % (Ohm)
R2 = 500;       % (Ohm)
R3 = 300;       % (Ohm)
L1 = 1;         % (H)
L2 = 2;         % (H)
C1 = 1e-3;      % (F)
C2 = 2e-3;      % (F)

I  = 5;         % (A)

%%%%%%%%%%%%%%%%
% Solve
%%%%%%%%%%%%%%%%
s = tf('s');                % Laplace operator

% Impedances
zr1   = R1;
zr2   = R2;
zr3   = R3;

zl1   = s*L1;
zl2   = s*L2;

zc1   = 1/(s*C1);
zc2   = 1/(s*C2);

zr2l2 = R2 + s*L2;

% Y matrix
Y = [ 1/zr1 + 1/zr3 + 1/zl1 + 1/zc1  -1/zr3
     -1/zr3                           1/zr3 + 1/zc2 + 1/zr2l2];

% I vector
I = [0 I/s]';

% Solve for voltage vector
V = inv(Y) * I;

% Display the 2 transfer functions that result
% Could have just omitted the ";" from the previous command
V

% Compute voltages
[v t] = impulse(V);

%%%%%%%%%%%%%%%%
% Plot
%%%%%%%%%%%%%%%%
% Plot both voltages together
figure(1); clf; hold on; grid on;
plot(t, v(:,1), 'k-', 'Linewidth', 3);
plot(t, v(:,2), 'b:', 'Linewidth', 3);
title('Node Voltages');
ylabel('Voltage (V)');
legend('Node 1', 'Node 2', 'Location', 'East');
set(gca, 'FontSize', 14);

% Re-Plot voltages at node 1
figure(2); clf; hold on; grid on;
plot(t, v(:,1), 'k-', 'Linewidth', 3);
title('Node  1 Voltage');
ylabel('Voltage (V)');
set(gca, 'FontSize', 14);

