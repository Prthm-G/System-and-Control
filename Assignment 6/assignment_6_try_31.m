% Changed to multiplying variables
% Assignment 6
% Eve Sankar
% 14425508
% Try 3

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

%Variables
CF = 10*F;
G_DCG = 300/G;
PA = -A;
PB = -3*B;
ZC = -5*C;
PD1 = 2*D*(-1+i);
PD2 = 2*D*(-1-i);
PE = -25*E;

% Q1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create a transfer function from the figure 2

EMS = zpk(ZC, [PA, PB, PD1, PD2], 1);
SEN = zpk([], PE, 1);
kf = 1/dcgain(SEN);
k_EMS = 1/dcgain(EMS);
dh = CF/(s+CF);

Q1.G = EMS*300/G*k_EMS;
Q1.GH = EMS*300/G*SEN*dh*kf*k_EMS;
Q1.Ku = margin(Q1.GH);
Q1.K = Q1.Ku/2;
Q1.X = Q1.K*Q1.G/(1+Q1.K*Q1.GH); 

Q1.Ess = dcgain(1/(1+Q1.K*Q1.GH))*100;

% Q2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% same transfer functions

Q2.G = Q1.G/s;
Q2.GH = Q1.GH/s;
Q2.Ku = margin(Q2.GH);
Q2.K = Q2.Ku/2;
Q2.X = Q2.K*Q2.G/(1+Q2.K*Q2.GH);
Q2.Ess = dcgain(1/(1+Q2.K*Q2.GH))*100;

% Q3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We want a dcgain(1/(1+Q3.K1*Q2.GH))*100 = 30
% s = 0

Q3.K1 = (1/0.3 - 1)/dcgain(Q1.GH)

% % We want 0% overshoot
% % Ess is min when K is max therefore K = Ku

Q3.K2 = Q1.K/7.35;
Q3.X2 = Q3.K2*Q1.GH/(1+Q3.K2*Q1.GH);
Q3.Ess_2 = dcgain(1/(1+Q3.K2*Q1.GH))*100

% done

% COW, Ess = 15ish%

% Q4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 10% OS

Q4.K1 = 0.22235*Q2.Ku;
Q4.X1 = Q4.K1*Q2.GH/(1+Q4.K1*Q2.GH);

% Q4.K1 = 21.05ish percent of Q2.Ku

% 0% overshoot
Q4.K = Q1.Ku/9.072;
Q4.X2 = Q4.K2*Q2.GH/(1+Q4.K2*Q2.GH);
figure(2)
step(Q4.X2)
stepinfo(Q4.X2)

if 0
    SN = 14425508
    makeMat341
end

% 
% % Eve prattles on
% % maybe (s+7*F) = (s+25*E). One of the poles is equal to 2CF
% % feedback path is composed of Dh and Kh
% %Q1.DH = 2*CF/(s+2*CF)
% %Q1.KH = 1/(dcgain(Hs))
% 
