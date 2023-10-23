% Note: All variables are cleared when this is run
clear all; clc;

%set student number variables
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

F0 = 100;
I0 = F0;

%%%%%%%%%%
% Q1
%%%%%%%%%%

%define system variables
syms X0 X1 X2 X3 s

%solve system of equations (output is a symbolic expression???)
Q1kcl0 = I0/s -X0/L0*s -X0 * s * C0 == (X0-X2)/R20 + (X0-X2)/(s*L20);
Q1kcl1 = (X2-X1)/(R21) - (X1)/(s*L1) -(X1 * s * C1) == (X1-X3)/R31;
Q1kcl2 = (0-X2*s*C2) + (X0-X2)/(R20) + (X0-X2)/(s*L20) == (X2-X3)/(s*L32) + (X2-X1)/(R21);
Q1kcl3 = (X1-X3)/(R31)+(X2-X3)/(s*L32) == X3*s*C3;
%solve for X0 X1 X2 X3
[v0_sym, v1_sym, v2_sym, v3_sym] = solve([Q1kcl0 Q1kcl1 Q1kcl2 Q1kcl3],[X0 X1 X2 X3]);

%convert symbolic expression to tf
ExpFun_v0 = str2func(regexprep(func2str(matlabFunction(v0_sym)), '\.([/^\\*])', '$1'));
v0 = tf(ExpFun_v0(tf('s')));

ExpFun_v1 = str2func(regexprep(func2str(matlabFunction(v1_sym)), '\.([/^\\*])', '$1'));
v1 = tf(ExpFun_v1(tf('s')));

ExpFun_v2 = str2func(regexprep(func2str(matlabFunction(v2_sym)), '\.([/^\\*])', '$1'));
v2 = tf(ExpFun_v2(tf('s')));

ExpFun_v3 = str2func(regexprep(func2str(matlabFunction(v3_sym)), '\.([/^\\*])', '$1'));
v3 = tf(ExpFun_v3(tf('s')));


%minreal and assign
s = tf('s');
v0 = minreal(v0);
v1 = minreal(v1);
v2 = minreal(v2);
v3 = minreal(v3);

%compute distsance  tranfer function
d0 = v0 * 1/s;
d1 = v1 * 1/s;
d2 = v2 * 1/s;
d3 = v3 * 1/s;

%compute tranfer function for d3
Q1tf_d1 = minreal(d1 / (I0/s))

figure(1);
impulse(Q1tf_d1);
grid on

Q1.G1 = Q1tf_d1;

Q1tf_d3 = tf(d3 / (I0/s))

figure(2);
impulse(Q1tf_d3);
grid on

Q1.G3 = Q1tf_d3;

%%%%%%%%%%
% Q2
%%%%%%%%%%

%seperating force  = current though K1 or distance between m0 and m2
%multiplied by spring constant 

Q2sfL1 = minreal((0-v1)/(s*L1));
Q2tf_sfL1 = minreal(Q2sfL1/(I0/s),0.0001)

figure(3);
impulse(Q2tf_sfL1);
grid on

Q2.G1 = Q2tf_sfL1;

Q2sfL32 = minreal((v2-v3)/(s*L32));
Q2tf_sfL32 = minreal(Q2sfL32/(I0/s),0.0001)

figure(4);
impulse(Q2tf_sfL32);
grid on

Q2.G32 = Q2tf_sfL32;


%%%%%%%%%%%
% Q4 and Q5
%%%%%%%%%%%

if 1
%set variables
Rw = numA/3;
Lw = (numB) * 10^(-3);

Jr = (numC / 10) * 10^(-3);
Br=(numD) * 10^(-3);

km = (50*numG)/1000;
Km = (50*numG)/1000;

Vin = 150; %input of 1 (for some reson i cant use imput 12/s then divide output by 12/s)

%define system variables (X -> Iw, Y -> W)
syms current omega s

%solve system of equations (output is a symbolic expression???)
Q45eq1 = current == (Vin/s - km*omega)/(Rw + s*Lw);
Q45eq2 = Km*current == (omega*Br + omega*Jr*s);
[Iw_sym, W_sym] = solve([Q45eq1 Q45eq2],[current omega]);

%convert symbolic expression to tf
ExpFun_W = str2func(regexprep(func2str(matlabFunction(W_sym)), '\.([/^\\*])', '$1'));
W = tf(ExpFun_W(tf('s')));
ExpFun_Iw = str2func(regexprep(func2str(matlabFunction(Iw_sym)), '\.([/^\\*])', '$1'));
Iw = tf(ExpFun_Iw(tf('s')));

s = tf('s');

%minreal and assign
tf_W = minreal(W / (Vin/s))
tf_Iw = minreal(Iw / (Vin/s))

figure(5);
step(tf_W);
grid on;

figure(6);
step(tf_Iw);
grid on;

Q3.G = tf_W;
Q4.G = tf_Iw;    
end

if 0
   SN = 90300815; 
end