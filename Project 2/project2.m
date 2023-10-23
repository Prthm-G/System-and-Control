%% Project 2
clear all; clc;
SN = 90300815;
A = 9 + 10;
B = 0 + 10;
C = 3 + 10;
D = 0 + 10;
E = 0 + 10;
F = 8 + 10;
G = 1 + 10;
H = 5 + 10;

% from Project 1
%p1DSPlot222(SN);
Peak = 13.6003;
FV = 10.995;
Kdc = FV;
Tr = 0.72 * 1E-3; % ms -> s
Tp = 1.52 * 1E-3; % ms -> s
Ts = 4.6 * 1E-3; % ms -> s
OS = (Peak - FV) / FV; % unitless

zeta = sqrt((log(OS) ^ 2) / (pi ^ 2 + log(OS) ^ 2));
beta = 1 - zeta ^ 2;
wnr = (pi - atan(beta / zeta)) / (Tr * beta);
wns = 4 / (zeta * Ts);
wn = (wnr + wns) / 2;

s = tf('s');
Gsp = Kdc * wn ^ 2 / (s ^ 2 + s * 2 * zeta * wn + wn ^ 2);

G1 = A * 7 / (s + B * 800);
G2 = C * 8 / (s + D * 700);
G3 = 10 ^ 5 / (s + E * 600);
G4 = F * 50 / (s + G * 500);
H1 = 4 / (s + H * 5);

% Set input output relationships and solve system of 8 equations
% Input U is 1
tfq3 = (G1 * G2 * G3 * G4 * H1 + G1 * G2 + G2 * G3) / (G2 * G3 * H1 + G3 * G4 * H1 + 1);
Khp = dcgain(tfq3) * 1E-3; % mV -> V
Dhp = tfq3 / dcgain(tfq3); % Since pure

% Actuator
Rw = A / 2; % ohm
Lw = B * 30 * 1E-6; % uH -> H
Km = C * 1E-3; % mNm/A -> Nm/A
Mm = (D + E) * 1E-3; % g -> kg
Bm = F / 30 * 1E-6; % uNms -> Nms
Jr = G / 15 * 1E-7; % gcm^2 -> kgm^2
Js = H / 5 * 1E-7; % gcm^2 -> kgm^2
Ms = A / 4 * 1E-3; % g -> kg
Bs = B / 3; % Ns/m
Na = 3 * 1E-2 / (2 * pi); % cm/turn -> m/rad
% Mechanism (3 fingers)
Jf = 3 * C / 3 * 1E-7; % gcm^2 -> kgm^2
Bf = 3 * D * 1E-3; % mNms -> Nms
Nf = 10 * (pi / 180) * 1E2; % deg/cm -> rad/m
Bt = 3 * E * 1E-3; % mNms -> Nms
Kt = 3 * F * 30 * 1E-3; % mNm -> Nm
Nt = 100 * 1E-3; % Nt = L6
Bl = 3 * G / 5; % Ns/m
Kl = 3 * H * 15; % N/m

Jeff = Jr + Js + Na^2 * (Mm + Ms) + Nf^2 * Na^2 * Jf;
Beff = Bm + Na^2 * Bs + Na^2 * Nf^2 * Bf;
wm_qm = (-Na^2 * Nf^2 * Nt^2 * Kl - Na^2 * Nf^2 * Kt) / (Nf^2 * Na^2 * Nt^2 * Bl + Na^2 * Nf^2 * Bt);
wm_qr = (Na^2 * Nf^2 * Kt) / (Nf^2 * Na^2 * Nt^2 * Bl + Na^2 * Nf^2 * Bt);
wm_wr = (Na^2 * Nf^2 * Bt) / (Nf^2 * Na^2 * Nt^2 * Bl + Na^2 * Nf^2 * Bt);

% x = [qr; wr; qm; Iw]
Amat = [0, 1, 0, 0;
         -Na^2 * Nf^2 * Kt / Jeff + Na^2 * Nf^2 * Bt / Jeff * wm_qr, -Beff / Jeff - Nf^2 * Na^2 * Bt / Jeff + Na^2 * Nf^2 * Bt / Jeff * wm_wr, Na^2 * Nf^2 * Kt / Jeff + Na^2 * Nf^2 * Bt / Jeff * wm_qm, Km / Jeff;
         wm_qr, wm_wr, wm_qm, 0;
         0, -Km / Lw, 0, -Rw / Lw];
Bmat = [0; 0; 0; 1 / Lw];

Cmat = [0, 0, 0, Km;
         Na, 0, 0, 0;
         (Na * Nf * 180) / pi, 0, 0, 0;
         Nf^2 * Na^2 * Nt^2 * Bl * wm_qr, Nf^2 * Na^2 * Nt^2 * Bl * wm_wr, Na^2 * Nf^2 * Nt^2 * Kl + Nf^2 * Na^2 * Nt^2 * Bl * wm_qm, 0;
         1, 0, 0, 0];

Dmat = [0;0;0;0;0];

output = ss(Amat, Bmat, Cmat, Dmat);

%% Q1: Feedback Path
Hsp = Dhp * 1E-3;
motorang = output(5) * 180 / pi;
GH = motorang * Hsp * Gsp;

CF = 15 * C;
Hs = Hsp / (Na * Nf);
Gs = Gsp * output(3);
Kh = 1 / dcgain(Hs);
Dh = CF / (s + CF); 
OL1 = Gs * Hs * Kh * Dh;
Ku1 = margin(OL1); 

% Set proportional (master) for P-control to 25% of Ku
K1 = Ku1 * 0.25;

CL1 = feedback(K1 * Gs, Hs * Kh * Dh);
info1 = stepinfo(CL1);

Q1.Kh = Kh;
Q1.Ndt = 1;
Q1.GOS = (info1.Peak - 1) * 100;
Q1.Ess = (1 - dcgain(CL1)) * 100;

%% Q2: FDD WS Filter
Nd = B;
n = 1 : Nd;
Wp = exp(-4 * (n - 1) / (Nd - 1));
W = Wp / sum(Wp);
Nhat = sum((n - 1) .* W);

Q2.W = W;
Q2.Ndhat = Nhat;

%% Q3: Controller Metrics
Ndt = Nhat + 0.5; % includes reading and writing
% Partial dynamics block: all 0s=inf, gains=unity
Dp = (1 / s) * (CF / (Ndt * s + CF));
p = pole(Dp)'; 
% Ultimate gain of control system, K=1, 0s=inf
Ku2 = margin(Dp * Gs * Kh * Hs * Dh);
% K=Ku2 for crossover f
[gm3, pm3, wxo3] = margin(Ku2 * Dp * Gs * Kh * Hs * Dh);

Q3.p = p;
Q3.Kappa = Ku2;
Q3.wxo = wxo3;

%% Q4: Initial Zeros
% To cancel 2 0s, one integral and one derivative controlled
% We can optimize center distance on wxo
% Then optimize separation distance
%{
[DzDp, zeros, maxPm] = maxPmRealZeros(wxo3, Ku2, Dp, Gs * Hs * Kh * Dh);
if (zeros(1) == zeros(2))
    [DzDp, zeros, maxPm] = maxPmComplexZeros(wxo3, Ku2, Dp, Gs * Hs * Kh * Dh);
end
%}
zeros = [(-16.1710 +27.1300i) (-16.1710 -27.1300i)];
maxPm =   108.4182;
DzDp = (0.1955*s^2 + 6.322*s + 195) / (2.168 * s^2 + 195 * s);

Ku4 = margin(DzDp * Gs * Hs * Kh * Dh);
K4 =  23.1320; % start at Ku/4
CL4 = feedback(K4 * DzDp * Gs, Hs * Kh * Dh);
[Gm4, Pm4, wxo4] = margin(K4 * DzDp * Gs * Kh * Hs * Dh);
% {
while (Pm4 < 30)
    K4 = K4 - 1e-6
    [Gm4, Pm4, wxo4] = margin(K4 * DzDp * Gs * Kh * Hs * Dh);
    display(Pm4);
end
%}
Q4.z = zeros;
Q4.PM = maxPm;
Q4.D = DzDp;
Q4.K = 14.1353;

%% Q5: Initial PID Gains

