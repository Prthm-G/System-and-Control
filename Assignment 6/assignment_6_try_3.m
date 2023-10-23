
% Assignment 6
clear all; clc;

s = tf('s');

SN = 90300815;
% Student Number
A = 19;
B = 10;
C = 13;
D = 10;
E = 10;
F = 18;
G = 11;
H = 15;

if 1
%Variables
G_DCH = 5/F;
G_DCG = 3/G;
PA = -A;
PB = -3*B;
ZC = -5*C;
PD1 = 2*D*(-1+i);
PD2 = 2*D*(-1-i);
PE = -15*E;

% Q1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create a transfer function from the figure 2

EMS = zpk(ZC, [PA, PB, PD1, PD2], 1);
SEN = zpk([], PE, 1);
kf = 1/dcgain(SEN);
k_EMS = 1/dcgain(EMS);
%dh = CF/(s+CF);

Q1.G = EMS*3/G*k_EMS;
Q1.GH = EMS*3/G*SEN*G_DCH*kf*k_EMS;
Q1.Ku = margin(Q1.GH);
Q1.Kh = 1/G_DCH;

figure(1)
step(Q1.GH);


%%%%%% Question 2 %%%%%%%%%
Q2.K = Q1.Ku/2;
Q2.Gcl = Q2.K*Q1.G/(1+Q2.K*Q1.GH);
Q2.Ess = dcgain(1/(1+Q2.K*Q1.GH))*100
Q2.Ts = 0.8170;
Q2.OS = 48.5564;
Q2.GOS = 8.96216426478;

figure(10)
step(Q2.Gcl);
stepinfo(Q2.Gcl)

% Q3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We want a dcgain(1/(1+Q3.K1*Q2.GH))*100 = 30
% s = 0
Q3.K = Q1.Ku/1.7263;
Q3.Gcl = Q3.K*Q1.G/(1+Q3.K*Q1.GH);
Q3.Ess = dcgain(1/(1+Q3.K*Q1.GH))*100;
Q3.Ts = 0.8170;

figure(20)
step(Q3.Gcl);
stepinfo(Q3.Gcl)




% Q4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Q4.K = Q1.Ku/1.00016;
Q4.Gcl = Q4.K*Q1.G/(1+Q4.K*Q1.GH);
Q4.Ess = dcgain(1/(1+Q4.K*Q1.GH))*100;
Q4.OS = 85.9199;
Q4.GOS = 67.059;

figure(20)
step(Q4.Gcl);
stepinfo(Q4.Gcl)


end

if 0
    SN=90300815;
end


