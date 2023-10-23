clear all; clc;
SN = 90300815; 

%%xfCheck222
% Student Number
A = 9 + 10;
B = 0 + 10;
C = 3 + 10;
D = 0 + 10;
E = 0 + 10;
F = 8 + 10;
G = 1 + 10;
H = 5 + 10;
s=tf('s');

CF = C*20;  % Hz
Ns = D;     % control cycles
Nd = D + E; % control cycles

% Table 2: Sensor
s = tf('s');
K = 100000;
G1 = s + (A*25);
H1 = (s/15) + B;
H2 = s^2 +(C*50*s) + (D*E*5000);

% Table 3: Motor
Rw = A/2;           % ohms
Lw = B*10*1e-3;     % mH -> H
Km = (C+D)*1e-2;    % Nm/A
Bm = (E+F)*1e-6;    % Nms/rad
Jm = G*3*1e-6;      % Nms^2/rad

% Table 4: Mechanism
Bg = D;         % Nms/rad
Bp = E*20;      % Nms/rad
J1 = F*2*1e-6;  % Nms^2/rad
J2 = G*2*1e-5;  % Nms^2/rad
J3 = H*2*1e-4;  % Nms^2/rad

%% Q1

%%xfDS222(SN);

finalValue = 1.86667 ;
peakValue = 2.30513;
Q1.Tp = 0.018; %seconds

kdc = finalValue;
os = (peakValue - finalValue)/finalValue;

[zeta,wn,xfer] = O2approx_Tp(kdc,os,Q1.Tp);

Q1.Ga = xfer;
%%step(Q1.Ga);

%% Q2

Q2.Hs = (G1* K) / ((1 + H1) * (1 + H2));

%% Q3
N = 3; %Gear Ratiotf
Q3.B = Bm + Bg/(N^2) + (Bg+Bp)/(N^4);
Q3.J = Jm + J1 + J2/(N^2) + J3/(N^4);
Q3.K = 0;

%% Q4

Q4.Ye = 1/(s*Lw + Rw);
Q4.Ym = 1/((Q3.J)*s + (Q3.B) + (Q3.K)/s);

Q4.Gm = feedback(Q4.Ye*Km*Q4.Ym, Km);
Q4.Gs = Q4.Gm*Q1.Ga;

%% Q5

gain = 1/dcgain(Q2.Hs);
Q5.GH = Q4.Gm * Q1.Ga * Q2.Hs * 1/N;

%% Q6
n = 1 : Nd;
Wp = exp(-4 * (n - 1) / (Nd - 1));
Wd = Wp / sum(Wp);
Q6.Ndhat = sum((n - 1) .* Wd); % 7.430471184884621

n = 1 : Ns;
Wp = exp(-4 * (n - 1) / (Ns - 1));
Ws = Wp / sum(Wp);
Q6.Nshat = sum((n - 1) .* Ws); % 3.735701392077635

Q6.Nhat = Q6.Ndhat+Q6.Nshat;

Q6.Kf = 1 / dcgain(Q2.Hs / N);
Q6.Df = CF / ((Q6.Nhat + 1) * s + CF);

%% Q7

Q7.Dp = getDpPID(-2*CF, Q6.Ndhat);  

oltf = Q5.GH*Q6.Kf*Q6.Df;

Q7.Kref = margin(oltf*Q7.Dp);
Q7.K0=Q7.Kref;

[GM, PM, wxo] = margin(Q7.Kref*oltf*Q7.Dp);

Q7.wxo = wxo;

%% Q8

%[Dz8, zreal, maxPM] = maxPmRealZeros(Q7.Wxo,Q7.Kref,Q7.Dp,oltf);
%[Dz8, zreal, maxPM] = maxPmComplexZeros(Q7.Wxo,Q7.Kref,Q7.Dp,oltf);

Dz8 = (0.1221*s^2 + 4.396*s + 131.1)/(s^2 + 131.1*s);
zreal = [-18.0026+27.3800i  -18.0026-27.3800i];
maxPM =  92.6436;

Q8.D=Dz8;
Q8.Z = zreal;
Q8.PM = maxPM;
%Q8.K = getPhaseMarginK(Q7.Kref,Dz8,oltf,5);

%% Q9
% get initial values for Ki,Kp and Kd

Kp = getKpPID(Q8.Z(1), Q8.Z(2), -2*CF, 1);
Kd = getKdPID(Q8.Z(1), Q8.Z(2), -2*CF, 1, Kp);
Ki = 1;
%intial values:
% Kp = 3.4463
% Kd = 6.0369
% Ki = 1

masterK = 0.8;
Kp_test = 5.9;
Kd_test = 3.4;
Ki_test = 0.9;
% get dynamics so we can plot and check if we meet RCGs

dynamicsq9 = getDynamicsPID(Kp_test,Ki_test,Kd_test,CF,1);
tfn9 = feedback(masterK* dynamicsq9 * Q1.Ga * Q4.Gm * 1/s, Q2.Hs * Q6.Kf * Q6.Df);
step(tfn9);
stepinfo(tfn9)




Q9.K = masterK;
Q9.Ki = Ki_test;
Q9.Kd = Kd_test;
Q9.Kp = Kp_test;
