% Assignment 7 Try 2
% Eve Sankar

clear all; clc;

s = tf('s');

% Student Number
A = 11;
B = 14;
C = 14;
D = 12;
E = 15;
F = 15;
G = 10;
H = 18;

% Variables
CF = 10*F;
G_DCG = 300/G;
PA = -A;
PB = -3*B;
ZC = -5*C;
PD1 = 2*D*(-1+i);
PD2 = 2*D*(-1-i);
PE = -25*E;
P = -2*CF;

%% from Asn 6

EMS = zpk(ZC, [PA, PB, PD1, PD2], 1);
SEN = zpk([], PE, 1);
kf = 1/dcgain(SEN);
k_EMS = 1/dcgain(EMS);
dh = CF/(s+CF);

Q1G = EMS*300/G*k_EMS;
Q1GH = EMS*300/G*SEN*dh*kf*k_EMS;

%% Q1

Z = PA;
Q1.Kd = (Z-P)/(P*Z); % from Katie's notes
Q1.D = 1 + Q1.Kd*(-P*s/(s-P));

%% Q2

Q2.DGH = Q1.D*Q1GH;
Q2.Ku = margin(Q2.DGH);
Q2.Ess = dcgain(1/(1+Q2.Ku*0.5*Q2.DGH))*100;
Q2.Ess99 = dcgain(1/(1+Q2.Ku*0.99*Q2.DGH))*100;
% makes sense, because steady state error 
% goes minimum the closer you get to Ku

%% Q3
% Need peak value < 1.2
% Ess < 32.7731%

Q3.Z = -33; % smallest Ess with this %33
Q3.D = 1+(Q3.Z-P)/(P*Q3.Z)*(-P*s/(s-P));
Q3.DGH = Q3.D*Q1GH;
Q3.Ku = margin(Q3.DGH);
Q3.K = Q3.Ku/1.82 % then we refine this %1.82
Q3.DG = Q3.D * Q1G;
Q3.X = Q3.K*Q3.DG/(1+Q3.K*Q3.DGH);
% then you go back to q3.z, but it was fine and the smallest

Q3.Ess = dcgain(1/(1+Q3.K*Q3.DGH))*100

%% Q4

Q4.DGH = Q1.D*Q1GH*1/s;
Q4.DG = Q1.D*Q1G*1/s;
Q4.Ku = margin(Q4.DGH);
Q4.K = Q4.Ku/2;
Q4.X = Q4.K*Q4.DG/(1+Q4.K*Q4.DGH);

%% Q5
% need overshoot less than 20%
% settling time as small as possible

Q5.Z = -9; % smallest ts with this
Q5.D = 1+(Q5.Z-P)/(P*Q5.Z)*(-P*s/(s-P));
Q5.DGH = Q5.D*Q1GH*1/s;
Q5.Ku = margin(Q5.DGH);
Q5.K = Q5.Ku/2.528 % then we refine this for overshoot
Q5.DG = Q5.D * Q1G*1/s;
Q5.X = Q5.K*Q5.DG/(1+Q5.K*Q5.DGH);
% then you go back to q3.z, but it was fine and the smallest

stepinfo(Q5.X)

if 0
    SN = 14425508;
    makeMat341
end

%% Neanderthal value checking
% at PA, ess = 28.1835
% K = q3.ku/1.61
% 
% Pa -100
% Ess = 23.8743
% 150
% 25
% 
% 33 19.9165
% 34 19.9191
% 32 19.9259
% 35 19.9323
% 
% now with q3.ku/1.82
% -31, ess = 21.9789
% -33, ess = 21.9442
% -35 ess 21.9612
% 
% q5
% 10 0.2663
% 16 0.5449
% 
% with q5.k = q5.ku/2.528
% -10, 0.4826
% 11 .5013
% 9 .4696

