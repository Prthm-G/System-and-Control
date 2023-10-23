%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solution to Assignment 7
% 2022/11/10
% Katie Seifert
% Student #68469311
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clean up workspace
clear all; clc;

% SN variable must contain Student Number
% This must be right or solution will not be graded
SN    = 68469311;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Variables

A = 6 + 10;
B = 8 + 10;
C = 4 + 10;
D = 6 + 10;
E = 9 + 10;
F = 3 + 10;
G = 1 + 10;
H = 1 + 10;

s = tf('s');
%Variables
CF = 10*F;
G_DCG = 300/G;
PA = -A;
PB = -3*B;
ZC = -5*C;
PD1 = 2*D*(-1+i);
PD2 = 2*D*(-1-i);
PE = -25*E;

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

Q1G = EMS*300/G*k_EMS;
Q1GH = EMS*300/G*SEN*dh*kf*k_EMS;
%Q1.Ku = margin(Q1.GH);
%Q1.K = Q1.Ku/2;
%Q1.X = Q1.K*Q1.G/(1+Q1.K*Q1.GH); 
%Q1.Ess = dcgain(1/(1+Q1.K*Q1.GH))*100;

%Q2
% same transfer functions

Q2G = Q1G/s;
Q2GH = Q1GH/s;
%pzmap(Q2GH)
%Q2Ku = margin(Q2GH);
%Q2K = Q2Ku/2;
%Q2X = Q2K*Q2G/(1+Q2K*Q2GH);
%Q2Ess = dcgain(1/(1+Q2K*Q2GH))*100;

%% Q1
%pzmap(Q1GH)
Z = PA;
Kd = (Z-P)/(P*Z)
Dd = 1+Kd*(-P*s/(s-P))
%for reference in case D and K arent seperate D = 1+(Z-P)/(PZ)*(-Ps/(S-P))
%pzmap(Q1.GH*(1+K*D))
Q1.Kd = Kd;
Q1.D = Dd;

%% Q2
Ku = margin(Dd*Q1GH)
K1 = Ku*.5;
K2 = Ku*.99;
Ess50 = dcgain(1/(1+K1*Dd*Q1GH))*100
Ess99 = dcgain(1/(1+K2*Dd*Q1GH))*100
Q2.Ku = Ku;
Q2.Ess = Ess50; %might need to be named Ess50
Q2.Ess99 = Ess99;

%% Q3
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
%{
%outdate loop
%loop to test things
t = 0;
pt = 10;
min = 199;
%got poles on line @ kt = 0.206, Z=17.3
while t < Ku
Kt = t;
Z = PA;
Kd = (Z-P)/(P*Z);
Dd = 1+Kd*(-P*s/(s-P));
X = Dd*Kt*Q1G/(1+Dd*Kt*Q1GH);
%pzmap(X)
%step(X)
Esst = dcgain(1/(1+Kt*Dd*Q1GH))*100;
S = stepinfo(X);
if ( S.Peak < 1.2 && Esst<Ess50)
    disp("\n")
    disp(S.Peak)
    disp(Kt)
    disp(Ku/t)
    disp(Z)
    disp(Esst)
    if (min > Esst)
        min = Esst;
    end
end
t = t+0.01*Ku;
end
disp(min)

%figure, hold on
%{
while t < Ku
Kt = t
Z = PA
Kd = (Z-P)/(P*Z)
Dd = 1+Kd*(-P*s/(s-P))
X = Dd*Kt*Q1G/(1+Dd*Kt*Q1GH)
step(X)
t = t + .1*Ku;
end
%}
%}

Kt = 0.1142
Z = -49.6000
Kd = (Z-P)/(P*Z);
Dd = 1+Kd*(-P*s/(s-P));
X = Dd*Kt*Q1G/(1+Dd*Kt*Q1GH)
Q3.K = Kt;
Q3.Z = Z;
Q3.X = X;

%% Q4
Z = PA;
Kd = (Z-P)/(P*Z)
Dd = 1+Kd*(-P*s/(s-P))
Ku = margin(Dd*Q2GH)
K50 = Ku/2
X = Dd*K50*Q2G/(1+Dd*K50*Q2GH)
Q4.Ku = Ku;
Q4.X = X;

%% Q5
%Requirements
% Overshoot < 20%
%Goal: Minimize Ts
%{
%loop to test things
t = 0;
pt = -10*16;
min = 1000;
kf = 0
zf = 0
while t < Ku
    Kt = t;
    while pt < -1
    %Kt = 0.0114; best result
    Z = pt;
    Kd = (Z-P)/(P*Z);
    Dd = 1+Kd*(-P*s/(s-P));
    X = Dd*Kt*Q2G/(1+Dd*Kt*Q2GH);
    %pzmap(X)
    %step(X)
    S = stepinfo(X);
    if ( S.Overshoot < 20)
        disp("/n")
        disp(S.Overshoot)
        disp(Kt)
        disp(Z)
        disp(S.SettlingTime)
        if (min > S.SettlingTime)
            min = S.SettlingTime;
            kf = Kt;
            zf = Z;
            disp("min=")
            disp(min)
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
Q5.K = 0.3502;
Q5.Z = -14.4;
Kd = (Q5.Z-P)/(P*Q5.Z)
Dd = 1+Kd*(-P*s/(s-P));
X = Dd*Q5.K*Q2G/(1+Dd*Q5.K*Q2GH)
Q5.X = X;
