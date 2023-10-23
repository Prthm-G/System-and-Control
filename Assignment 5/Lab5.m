% lab 5

clc;
clear;

SN = 90300815;
numA = 10 + 9;
numB = 10 + 0;
numC = 10 + 3;
numD = 10 + 0;
numE = 10 + 0;
numF = 10 + 8;
numG = 10 + 1;
numH = 10 + 5;

%Set variables
M0 = numA/5;
M1 = numB/10;
M2 = numC/10;
M3= numD/5;
C0 = M0;
C3 = M3;
C1 = M1;
C2 = M2;

B20 = numE/2;
B21 = numF/3;
B31 = numG/4;
R20 = 1/B20;
R21 = 1/B21;
R31 = 1/B31;

K0 = numA;
K1 = numB;
K20 = numC;
K32 = numD/3;
L0 = 1/K0;
L1 = 1/K1;
L20 = 1/K20;
L32 = 1/K32;

F0 = 1;
I0 = F0;


% Q1

% dx = [dd0, dd1, dd2, dd3, dv0, dv1, dv2, dv3];
% x = [d0, d1, d2, d3, v0, v1, v2, v3];
% u = [F0];

s = tf('s');

A(1, :) = [0, 0, 0, 0, 1, 0, 0, 0];
A(2, :) = [0, 0, 0, 0, 0, 1, 0, 0];
A(3, :) = [0, 0, 0, 0, 0, 0, 1, 0];
A(4, :) = [0, 0, 0, 0, 0, 0, 0, 1];
A(5, :) = [-K20/M0, 0, K20/M0, 0, -(B20/M0 + B0/M0), 0, B20/M0, 0];
A(6, :) = [0, -(K1/M1 + K21/M1), K21/M1, 0, 0, -B31/M1, 0, B31/M1];
A(7, :) = [K20/M2, K21/M2, -(K32/M2 + K20/M2 + K21/M2), K32/M2, B20/M2, 0, -B20/M2, 0];
A(8, :) = [0, 0, K32/M3, -K32/M3, 0, B31/M3, 0, -B31/M3];

B = zeros(8, 1);
B(5, 1) = 300/M0;

Q1.A = A;
Q1.B = B;

% Q2

% y = [d3, FK20];
% x = [d0, d1, d2, d3, v0, v1, v2, v3];
% u = [F0];

C = [0, 0, 0, 1, 0, 0, 0, 0; K20, 0, -K20, 0, 0, 0, 0, 0];

D = [0; 0];

Q2.C = C;
Q2.D = D;

% COW

[a, b] = ss2tf(A, B, C, D);

xfer1 = tf(a(1, :), b);
xfer2 = tf(a(2, :), b);

% Looks the same
%step(xfer1, xfer2)

% Q3

Rw = 18/3;
Lw = 16e-3;
K = 18/2;
Vin = 150;
M0 = 18+16;
R0 = 1/(11+10);
N = 2;

clear A;
A(1, :) = [0, 0, 0, 0, 1, 0, 0, 0, 0];
A(2, :) = [0, 0, 0, 0, 0, 1, 0, 0, 0];
A(3, :) = [0, 0, 0, 0, 0, 0, 1, 0, 0];
A(4, :) = [0, 0, 0, 0, 0, 0, 0, 1, 0];
A(5, :) = [-(K20/M0 + 1/(M0*R0)), 0, K20/M0, 0, -(B20/M0 + B0/M0), 0, B20/M0, 0, K/(2*M0)];
A(6, :) = [0, -(K1/M1 + K21/M1), K21/M1, 0, 0, -B31/M1, 0, B31/M1, 0];
A(7, :) = [K20/M2, K21/M2, -(K32/M2 + K20/M2 + K21/M2), K32/M2, B20/M2, 0, -B20/M2, 0, 0];
A(8, :) = [0, 0, K32/M3, -K32/M3, 0, B31/M3, 0, -B31/M3, 0];
A(9, :) = [0, 0, 0, 0, -K/Lw, 0, 0, 0, -Rw/Lw];

B = zeros(9, 1);
B(9, 1) = 1/Lw;

C = [0, 0, 0, 1, 0, 0, 0, 0, 0];

D = [0];

[num, den] = ss2tf(A, B, C, D);

xfer3 = tf(num, den);
xfer3 = minreal(xfer3);

%step(xfer3)

Q3.G = xfer3;

% Q4

B = zeros(9, 2);
B(9, 1) = 1/Lw;

C = [0, 0, 0, 1, 0, 0, 0, 0, 0];

D = [0, 0];

SS = ss(A, B, C, D);
xfer1 = tf(SS(1));

B = [0, 0; 0, 0; 0, 0; 0, 0; 0, -1; 0, -1; 0, -1; 0, -1; 0, 0];

SS = ss(A, B, C, D);
xfer2 = tf(SS(2));

%step(xfer2)

Q4.Gv = xfer1;
Q4.Gg = xfer2;





















