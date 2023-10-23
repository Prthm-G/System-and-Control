%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solution to Assignment 8
% 2022/11/18
% Katie Seifert
% Student #68469311
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clean up workspace
clear all; clc;

% SN variable must contain Student Number
% This must be right or solution will not be graded
SN    = 90300815;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Variables

A = 9 + 10;
B = 0 + 10;
C = 3 + 10;
D = 0 + 10;
E = 0 + 10;
F = 8 + 10;
G = 1 + 10;
H = 5 + 10;

s = tf('s');
%Variables
CF = 10*G;
G_DCG = 3/G;
G_DCH = 5/F;
PA = -A;
PB = -3*B;
ZC = -5*C;
PD1 = 2*D*(-1+i);
PD2 = 2*D*(-1-i);
PE = -15*E;

%new derivates
P = -2*CF;


%% Asn6
%Q1
% create a transfer function from the figure 2

EMS = zpk(ZC, [PA, PB, PD1, PD2], 1);
SEN = zpk([], PE, 1);
kf = 1/dcgain(SEN);
k_EMS = 1/dcgain(EMS);
dh = CF/(s+CF);

Q1G = EMS*3/G*k_EMS;
Q1GH = EMS*3/G*SEN*dh*kf*k_EMS;
%pzmap(Q1GH)
%Q1.Ku = margin(Q1.GH);
%Q1.K = Q1.Ku/2;
%Q1.X = Q1.K*Q1.G/(1+Q1.K*Q1.GH); 
%Q1.Ess = dcgain(1/(1+Q1.K*Q1.GH))*100;

%Q2
% same transfer functions

%Q2G = Q1G/s;
%Q2GH = Q1GH/s;
%pzmap(Q2GH)
%Q2Ku = margin(Q2GH);
%Q2K = Q2Ku/2;
%Q2X = Q2K*Q2G/(1+Q2K*Q2GH);
%Q2Ess = dcgain(1/(1+Q2K*Q2GH))*100;

%% Q1
%pzmap(Q1GH)
%Requirements
% Peak value < 1.2
% Ess < Q2.Ess50
%Goal: Minimize Ess
%{
 Newer loop (tests Z and K)
%loop to test things
t = 0;
pt = -10*16;
min = 1000;
kf = 0;
zf = 0;
while t < Ku
    Kt = t;
    while pt < -1
    Z = pt;
    Kd = (Z-P)/(P*Z);
    Dd = 1+Kd*(-P*s/(s-P));
    X = Dd*Kt*Q1G/(1+Dd*Kt*Q1GH);
    %pzmap(X)
    %step(X)
    Esst = dcgain(1/(1+Kt*Dd*Q1GH))*100;
    S = stepinfo(X);
    if ( S.Peak < 1.2 && Esst<Ess50)
        if (min > Esst)
            min = Esst;
            kf = Kt;
            zf = Z;
        end
    end
    pt = pt + 0.1*16;
    end
    pt = -10*16;
    t = t+0.01*Ku;
end
disp("final values")
disp(min)
disp(kf)
disp(zf)
%}
Z = PA;
Kp = 1/(-Z);
Kd = (Z-P)/(P*Z);
Dd = Kd*(-P*s/(s-P));
Q1.Kp = Kp;
Q1.D = Dd;
Q1.K = margin(Dd*Q1GH);
Q1.X = Dd*Q1G/(1+Dd*Q1GH);
stepinfo(Dd*Q1GH)
[Gm,Pm,Wcg,Wcp] = margin(Dd*Q1GH)

figure(1)
step(Dd*Q1GH);

