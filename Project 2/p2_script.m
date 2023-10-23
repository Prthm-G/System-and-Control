%% Initial Setup 
clear; clc;
SN=90300815;
A=9+10;
B=0+10;
C=3+10;
D=0+10;
E=0+10;
F=8+10;
G=1+10;
H=5+10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%                    Project 1                      %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Variables: make everything SI units (g ->kg, cm->m, deg->rad)
Rw = A/2; % ohms
Lw = B*30e-6; % H
Km = C*1e-3; % Nm/A

Mm = (D + E)*1e-3; % Kg
Bm = (F/30)*1e-6; % Nms
Jr = (G/15)*1e-7; % Kg-m^2

Js = (H/5)*1e-7; % Kg-m^2
Ms = (A/4)*1e-3; % Kg
Bs = B/3; % Ns/m

Na = 3*1e-2/(2*pi); % m/rad

% Multiply finger and tine parameters by 3 since there are three fingerss
Jf = 3*(C/3)*1e-7; % Kg-m^2
Bf = 3*D*1e-3; % Nms
Nf = 10*pi/(180)*1e2; % rad/m

Bt = 3*E*1e-3; % Nms
Kt = 3*F*30e-3; % Nm

BI = 3*G/5; % Ns/m
KI = 3*H*15; % N/m

L6 = 100e-3; % m
Nt = L6; % m/rad

%% Q1: Black-Box Specs
Tp = 1.52; % ms
Tr = 0.72; % ms
Ts = 4.6; % ms

kDC = 10.9995; % DC gain and final value
OS = (13.6003-10.995)/10.995 ; % Overshoot
pOS = OS*100; % Percent overshoot

%% Q2: Second Order Approx
% Get natural frequencies for Tr and Ts, then take the mean to get wn
s = tf('s');
zeta = sqrt((log(OS))^2/(pi^2 + (log(OS))^2));
beta = sqrt(1-zeta^2);

wn_Tr = 1/(Tr*beta*1e-3)*(pi-atan(beta/zeta));
wn_Ts = 4/(zeta*Ts*1e-3);
wn = (wn_Tr + wn_Ts)/2; % average natural frequency

Q2G = kDC*(wn^2/(s^2+2*zeta*wn*s+wn^2)); % Second order approximation
t = 0:5/1000:5; % something wrong with how this is defined...  

%% Q3: Analog Position Sensor (BD manipulation)
% G1 = (A*7)/(s+B*800);
% G2 = (C*8)/(s+D*700);
% G3 = (1e5)/(s+E*600);
% G4 = (F*50)/(s+G*500);
% 
% H1 = 4/(s+H*5);
% 
% load sensor_ss; % State space generated from linear model of BD in Simulink
% D3 = tf(LinearAnalysisToolProject.LocalVariables.Value); % Simulink transfer function
% Kdc = dcgain(D3)/1000; % Getting DC gain from Simulink transfer function (convert mV to V)
% 
% 
% D3_manip = (G1*G2+G3*G2)/(1+H1*G3*G2); % derived from hand using BD manipulation
%% %% Question 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Weifeng's code
G1=A*7/(s+B*800);
G2=C*8/(s+D*700);
G3=10^(5)/(s+E*600);
G4=F*50/(s+G*500);
H1=4/(s+H*5);

%this is in mv/degree 
Q3_TF=((G1*G2*G3*G4*H1)+(G1*G2)+(G2*G3))/((G2*G3*H1)+(G3*G4*H1)+1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s = tf('s');


G3 = Q3_TF;
num = ((G1*G2*G3*G4*H1)+(G1*G2)+(G2*G3));
den = ((G2*G3*H1)+(G3*G4*H1)+1);
tfn3 = tf(num,den);

kdc3 = dcgain(tfn3);
Q3.Kdc = kdc3*10^(-3);

D3 = tfn3/kdc3;
Q3.D = D3;

figure(30);
step(Q3.D);
%% Q4-7: Mech/Electro-mech analysis and transfer functions
J_eff = Jr + Js + Na^2*(Mm + Ms) + Nf^2*Na^2*Jf;
B_eff = Bm + Na^2*Bs + Na^2*Nf^2*Bf;
wm_qm = (-Na^2*Nf^2*Nt^2*KI-Na^2*Nf^2*Kt)/(Nf^2*Na^2*Nt^2*BI + Na^2*Nf^2*Bt);
wm_qr = (Na^2*Nf^2*Kt)/(Nf^2*Na^2*Nt^2*BI + Na^2*Nf^2*Bt);
wm_wr = (Na^2*Nf^2*Bt)/(Nf^2*Na^2*Nt^2*BI + Na^2*Nf^2*Bt);

% x = [qr; wr; qm; Iw]
A_mat = [0, 1, 0, 0;
         -Na^2*Nf^2*Kt/J_eff + Na^2*Nf^2*Bt/J_eff*wm_qr, -B_eff/J_eff - Nf^2*Na^2*Bt/J_eff + Na^2*Nf^2*Bt/J_eff*wm_wr, Na^2*Nf^2*Kt/J_eff + Na^2*Nf^2*Bt/J_eff*wm_qm, Km/J_eff;
         wm_qr, wm_wr, wm_qm, 0;
         0, -Km/Lw, 0, -Rw/Lw];
B_mat = [0; 0; 0; 1/Lw];

C_mat = [0, 0, 0, Km;
         Na, 0, 0, 0;
         (Na*Nf*180)/pi, 0, 0, 0;
         Nf^2*Na^2*Nt^2*BI*wm_qr, Nf^2*Na^2*Nt^2*BI*wm_wr, Na^2*Nf^2*Nt^2*KI + Nf^2*Na^2*Nt^2*BI*wm_qm, 0;
         1, 0, 0, 0
        ];
D_mat = [0;0;0;0;0];

% C_mat = [1, 0, 0, 0;
%          0, 1, 0, 0;
%          0, 0, 1, 0;
%          0, 0, 0, 1
%          ];

sys = ss(A_mat, B_mat, C_mat, D_mat);
sys_tf = tf(sys);

%% Q8: Loop TF
% Input = voltage from amplifier
% Output = sensor voltage
motor_angle = sys_tf(5)*180/pi;
Hs = D3/1000;
GH = motor_angle * Hs * Q2G;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%                    Project 2                      %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Maybe we need to just use the full amplifier with motor angle transfer
% function as the system instead of just the finger angle transfer function


%% Q1: Feedback Path
CF = 15*C; % Control frequency
Hs = Hs/(Na*Nf); % We want to convert the finger angle to the motor angle  in the feedback path so
% we need to add this conversion in the sensor block

% The dimension we wish to control is finger angle qf, so this is what our
% "system" block will be (we also need to include the voltage amplifier for
% the motor since we need to actually power the motor)
Gs = Q2G*sys_tf(3); % Finger angle from Project 1 with the voltage amplifier needed for the motor voltage
% We have Hs = D3/1000 which was the sensor transfer function in V/deg
Kh = 1/dcgain(Hs); % we want our Kh to cancel the sensor gain
Dh = CF/(s+CF); % Accounting for read/write delays
% We do not need to add a weighted sum filter to the sensor, so the Ndt
% coefficient of the Dh feedback path delay will just be due to reading and
% writing

OpenLoop1 = Gs*Hs*Kh*Dh;
Ku1 = margin(OpenLoop1);

K1 = Ku1*0.25; % Set proportional (master, for P-control) to 25% of the ultimate gain Ku

ClosedLoop1 = feedback(K1*Gs,Hs*Kh*Dh);
stats1 = stepinfo(ClosedLoop1);
GOS1 = (stats1.Peak-1)*100;
Ess1 = (1 - dcgain(ClosedLoop1))*100;

Q1.Kh = Kh;
Q1.Ndt = 1;
Q1.GOS = GOS1;
Q1.Ess = Ess1;

% Now we make a PID controller to control the gripper

%% Q2: FDD WS Filter
% if we filter feedback path, we need Nhat + 1 (read/write)
% if we filter forward path, we need Nhat + 0.5 (derivative)
Nd = B;
n = 1:Nd;
Wp = exp(-4*(n-1)/(Nd-1));
W = Wp/sum(Wp);
Nhat = sum((n-1).*W);

Q2.W = W;
Q2.Ndhat = Nhat;
Ndt = Nhat + 0.5; % we do Ndt + 0.5 since this is a forwad path derivative weighted sum
% which does not ened to take into account reading and writin

% We will add this delay to the forward path derivative delay block to
% filter the finite difference derivative 

%% Q3: Controller Metrics
Dp = (1/s)*(CF/(Ndt*s + CF)); % partial dynamics of Dynamics block -> all zeros set to inf, gains set to unity
% We add the Ndt coefficient from the weighted sum on the derivative block
p = pole(Dp); % poles of the controller (using partial dynamics)
K0 = margin(Dp*Gs*Kh*Hs*Dh); % finding ultimate gain of control system with K=1, z=inf (partial dynamics)
[Gm1, Pm1, wxo1] = margin(K0*Dp*Gs*Kh*Hs*Dh); % now, set K=K0 and find the crossover frequency

Q3.p = p';
Q3.K0 = K0;
Q3.wxo = wxo1;

%% Q4: Initial Zeros
% Now, we cancel two zeros, one with the integral controller and one with
% the derivative controller
% Both zeros should be based off of wxo
% Place zeros about wxo, with some center distance
% Vary center distance until we find an optimal
% Then, vary seperation distance at the last 
%[dynamics, zeros, maxPm] = maxPmRealZeros(wxo1, K0, Dp, Gs*Hs*Kh*Dh);
%if(zeros(1) == zeros(2))
%    [dynamics, zeros, maxPm] = maxPmComplexZeros(wxo1, K0, Dp, Gs*Hs*Kh*Dh);
%end

%% 
%DzDp = dynamics; % dynamics is returned as the full PID dynamics based off
% the partial dynamics that was input to the function and the zeros that
% were returned

zeros = [(-16.0453 +26.9000i) (-16.0453 -26.9000i)];
maxPm =   108.2258;
DzDp = (0.1988*s^2 + 6.378*s + 195) / (2.168 * s^2 + 195 * s);
Q4.z = zeros;
Q4.PM = maxPm;
Q4.D = DzDp;

K04 = margin(DzDp*Gs*Hs*Kh*Dh); % find initial ultimate gain

% Now we tune our master gain to create a phase margin of 30deg
K4 = 23.1320;
ClosedLoop4 = feedback(K4*DzDp*Gs, Hs*Kh*Dh);
%%x1=masterGainTune(ClosedLoop4,30)
[Gm4, Pm4, wxo4] = margin(K4*DzDp*Gs*Kh*Hs*Dh);
% test = getK4PhaseMargin(K0,DzDp,Gs*Hs*Kh*Dh,30);
Q4.K = 23.1320;

%% Q5: Initial PID gains
% two poles, one at zero and one at CF/Ndt (use poles from controller)
Kp = 1/(p(2)) -(zeros(1)+zeros(2))/(zeros(1)*zeros(2));
Ki = 1;
Kd = 1/(zeros(1)*zeros(2)) + Kp/p(2);

% normalize the gains based off of K4 that made Pm=30deg so that we can
% start with a master gain of K=1 for the tuning
Q5.Kp = Kp*K4;
Q5.Ki = Ki*K4;
Q5.Kd = Kd*K4;

%% Q6: Meet RCGs (heuristic tuning)
% REQUIREMENTS
% GOS < 5%
% Ts < 100 ms
% Ess = 0%
%%%%%%%%%%%%%
% Goals
% GOS as small as possible
% Ts as small as possible

Kp6 = 0.75;
Ki6 = 16;
Kd6 = 0.0009;
K6 = 0.61;

D6 = Kp6 + Ki6/s + Kd6*(CF*s/(Ndt*s+CF));

ClosedLoop6 = feedback(K6*D6*Gs, Hs*Kh*Dh);
stats6 = stepinfo(ClosedLoop6);
GOS6 = (stats6.Peak - 1)*100; % Goal overshoot in percent
Ess6 = (dcgain(ClosedLoop6) - 1)*100; % Steady state error in percent
Ts6 = stats6.SettlingTime*1000; % Settling time in ms
step(ClosedLoop6);
xlim([0 0.3]);
Q6.K = K6;
Q6.Kp = Kp6;
Q6.Ki = Ki6;
Q6.Kd = Kd6;







