%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solution to Assignment 7
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
Nh = G;


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
Q1.wd = 1/CF;


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
Z = PA;
Kd = (Z-P)/(P*Z)
Dd = 1+Kd*(-P*s/(s-P))
%for reference in case D and K arent seperate D = 1+(Z-P)/(PZ)*(-Ps/(S-P))
%pzmap(Q1.GH*(1+K*D))no 
Q1.Kd = Kd;
Q1.D = Dd;

%% Q2
Ku = margin(Dd*Q1GH)
%K1 = Ku*.5;
%K2 = Ku*.99;
%Ess50 = dcgain(1/(1+K1*Dd*Q1GH))*100
%Ess99 = dcgain(1/(1+K2*Dd*Q1GH))*100
%Q2.Ku = Ku;
%Q2.Ess = Ess50; %might need to be named Ess50
%Q2.Ess99 = Ess99;
Q1.Ku = margin(Q1GH);
Q2.K = Q1.Ku/2;
Q2.X = Q2.K*Q1G/(1+Q2.K*Q1GH); 
Q2.Ess = dcgain(1/(1+Q2.K*Q1GH))*100;
%% Q3

%Requirements
% Peak value =1
% Ess < Q2.Ess50
%Goal: Minimize Ess

%Newer loop (tests Z and K)
%loop to test things


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

%Kt = 0.1142
%Z = -49.6000
%Kd = (Z-P)/(P*Z);
%Dd = 1+Kd*(-P*s/(s-P));
%X = Dd*Kt*Q1G/(1+Dd*Kt*Q1GH)
%Q3.K = Kt;
%Q3.Z = Z;
%Q3.X = X;

Q3.Ndt = Nh+1;
Q3.wd = Q3.Ndt/CF;

%% Q4
Z = -40;
Kd = (Z-P)/(P*Z)
Dd = 1+Kd*(-P*s/(s-P))
Ku = margin(Dd*Q1GH)
%K50 = Ku/2
%X = Dd*K50*Q2G/(1+Dd*K50*Q2GH)
Q4.Ku = Ku;
%Q4.X = X;
Q4.K=6.8528;
Q4.Ess=34.8559;


Q5.K = 6.3292;
Q5.Z = PA;
Q5.Kd = (Q5.Z-P)/(P*Q5.Z)
Dd = 1+Kd*(-P*s/(s-P));
%X = Dd*Q5.K*Q2G/(1+Dd*Q5.K*Q2GH)
%Q5.X = X;
Q5.Ess=36.3292;

t = 0;
pt = -10*16;
min = 1000;
kf = 0;
zf = -20;
while t < Ku
    Kt = t;
    Z = -20;
    Kd = (Z-P)/(P*Z);
    Dd = 1+Kd*(-P*s/(s-P));
    X = Dd*Kt*Q1G/(1+Dd*Kt*Q1GH);
    %pzmap(X)
    %step(X)
    Esst = dcgain(1/(1+Kt*Dd*Q1GH))*100;
    S = stepinfo(X)
    if ( S.Peak < 1 && Esst<Q2.Ess)
        if (min > Esst)
            min = Esst;
            kf = Kt;
            zf = Z;
        end
    end
    t = t+0.01*Ku;
end
disp("final values")
disp(min)
disp(kf)
disp(zf)

Q6.K = 6.4262;
Q6.Z = PA;
Q6.Kd = (Q6.Z-P)/(P*Q6.Z)
Dd = 1+Kd*(-P*s/(s-P));
%X = Dd*Q5.K*Q2G/(1+Dd*Q5.K*Q2GH)
%Q5.X = X;
Q6.Ess=36.3292;







%% Q2
Ku = margin(Dd*Q1GH)
t = 0;
pt = -10*16;
min = 1000;
kf = 0
zf = 0
while t < Ku
    Kt = t;
    while pt < -1
    Z = pt;
    Kp = 1/(-Z);
    Dd = (s-Z)/(s*(-Z));
    X = Dd*Kt*Q1G/(1+Dd*Kt*Q1GH);
    S = stepinfo(X);
    if ( S.Overshoot <= 10)
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
%% Q2 Final answers
Kp = 1/(-zf);
Dd = (s-zf)/(s*(-zf));
X = Dd*kf*Q1G/(1+Dd*kf*Q1GH);
Q2.K = kf;
Q2.Z = zf;
Q2.X = minreal(Dd*kf*Q1G/(1+Dd*kf*Q1GH));

%% Q3
%D = -p/(z1*z2)*(s-z1)*(s-z2)/s(s-p)
%Kd = (Z-P)/(P*Z)
%Kp = 1/(-Z)

Z1 = PA
Z2 = PB

Dd = -P/(Z1*Z2)*(s-Z1)*(s-Z2)/(s*(s-P))
Kp = 1/P-(Z1+Z2)/(Z1*Z2)
Kd = 1/(Z1*Z2)+Kp/P
Q3.Kp = Kp;
Q3.Kd = Kd;
Q3.D = Dd;
%% Q4
Ku = margin(Dd*Q1GH)
K = 0.5*Ku
X = minreal(Dd*K*Q1G/(1+Dd*K*Q1GH))
Q4.Ku = Ku;
Q4.X = X;

%% Q5
Z1 = PD1
Z2 = PD2

Dd = -P/(Z1*Z2)*(s-Z1)*(s-Z2)/(s*(s-P));
Kp = 1/P-(Z1+Z2)/(Z1*Z2);
Kd = 1/(Z1*Z2)+Kp/P;
Ku = margin(Dd*Q1GH)
K = 0.5*Ku
X = minreal(Dd*K*Q1G/(1+Dd*K*Q1GH))
Q5.Ku = Ku;
Q5.X = X;

%% Q6
%{
min = 1000;
kf = 0;
z1f = 0;
z2f = 0;
for Kt=0:0.01*Ku:Ku
    for zt1 = 10*PA:-0.1*PA:-1
        for zt2 = 10*PB:-0.1*PB:-1
            Kp = 1/(-zt1);
            Kd = P/zt2;
            Dd = -P/(zt1*zt2)*(s-zt1)*(s-zt2)/(s*(s-P));
            X = Dd*Kt*Q1G/(1+Dd*Kt*Q1GH);
            S = stepinfo(X);
            if ( S.Overshoot <= 10)
                disp("/n")
                disp(S.Overshoot)
                disp(Kt)
                disp(Z)
                disp(S.SettlingTime)
                if (min > S.SettlingTime)
                    min = S.SettlingTime;
                    kf = Kt;
                    z1f = zt1;
                    z2f = zt2;
                    disp("min=")
                    disp(min)
                end
            end
        end
    end
end
%}
min = 1000;
kf = 0;
zf = 0;
for Kt=0:0.01*Ku:Ku
    Ztr = PA
    Zti = 0
    while Ztr+Zti > PD1
    %for Zt = 10*PA:-0.1*PA:-1
        Zt1 = Ztr+Zti
        Zt2 = Ztr-Zti
        Dd = -P/(Zt1*Zt2)*(s-Zt1)*(s-Zt2)/(s*(s-P));
        Kp = 1/P-(Zt1+Zt2)/(Zt1*Zt2);
        Kd = 1/(Zt1*Zt2)+Kp/P;
        X = Dd*Kt*Q1G/(1+Dd*Kt*Q1GH);
        S = stepinfo(X);
        if ( S.Overshoot <= 10)
            disp("/n")
            disp(S.Overshoot)
            disp(Kt)
            disp(Zt1)
            disp(Zt2)
            disp(S.SettlingTime)
            if (min > S.SettlingTime)
                min = S.SettlingTime;
                kf = Kt;
                zf1 = Zt1;
                zf2 = Zt2;
                disp("min=")
                disp(min)
            end
        end
        Ztr = Ztr + 0.1*PA;
        Zti = Zti + 0.1*i;
    end
end
disp("final values")
disp(min)
disp(kf)
disp(zf1)
disp(zf2)

%% Q6 Final answers
Dd = -P/(zf1*zf2)*(s-zf1)*(s-zf2)/(s*(s-P))
Kp = 1/P-(zf1+zf2)/(zf1*zf2)
Kd = 1/(zf1*zf2)+Kp/P
X = Dd*kf*Q1G/(1+Dd*kf*Q1GH)
zf1
zf2
Q6.K = kf;
Q6.Z = [zf1 zf2];
Q6.X = minreal(X);





