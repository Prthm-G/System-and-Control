% Elec 341 ASN3 Script
% Author: Brandon Just
clear; clc;
SN = createSN(10435667);

% Values given in assignment
[M0,M3] = deal(SN.A); [M1,M2] = deal(SN.B);
B0 = SN.C; B20 = SN.E; B31 = SN.F; B1 = SN.D;
[K1, K32] = deal(SN.G); [K20, K21] = deal(SN.H);
F0 = 300;

% Electrical Equivlances
I = F0;
C0 = M0; C1 = M1; C2 = M2; C3 = M3;
R0 = 1/B0; R20 = 1/B20; R31 = 1/B31; R1 = 1/B1;
L1 = 1/K1; L32 = 1/K32; L20 =1/K20; L21 = 1/K21; 

syms V0 V1 V2 V3 s
eq1 = V0/R0+V0*C0*s + V0/(s*L20) - V2/(s*L20)+V0/(R20) - V2/R20 == I;
eq2 = V2/(s*L20)-V0/(s*L20)+V2/R20-V0/R20+V2*s*C2 +V2/(s*L21) - V1/(s*L21) -V3/(s*L32) + V2/(s*L32) == 0;
eq3 = V1/(s*L1)+V1*s*C1 - V2/(s*L21) + V1/(s*L21) - V3/R31 + V1/R31 == 0;
eq4 = V3*s*C3 + V3/R31 - V1/R31 + V3/(s*L32) - V2/(s*L32) == 0;
[A,B] = equationsToMatrix([eq1,eq2,eq3,eq4],[V0, V1, V2, V3]);
Voltages = linsolve(A,B);
Q1.G = sym2tf((Voltages(4)/s)/300);

%% Q2
i20 = (Voltages(1)-Voltages(3))/(s*L20); % NOTE. CHECK CURRENT DIRECTION
Q2.G = sym2tf(i20/300);

%% Q3
eq3 = V1/(R1)+V1*s*C1 - V2/(s*L21) + V1/(s*L21) - V3/R31 + V1/R31 == 0;
[A,B] = equationsToMatrix([eq1,eq2,eq3,eq4],[V0, V1, V2, V3]);
Voltages = linsolve(A,B);
Q3.G = sym2tf((Voltages(4)/s)/300);

%% Q4
Rw = 1 + SN.A/10;
Lw = (100 + 10*SN.B)*10^(-6);
Jr = (SN.C/10)*10^(-6);
Br = (SN.D+SN.E+SN.F)*10^(-6);
Km = (10+SN.G)*10^(-3);
Kmp = (10+SN.G)*10^(-3);
Vin = 12;
Rr = 1/Br;

syms w Iw
eq1 = Kmp*Iw == w/(1/(Jr*s))+w/Rr;
eq2 = (Vin-Km*w)/(s*Lw+Rw) == Iw;
[A,B] = equationsToMatrix([eq1,eq2],[w, Iw]);
Currents = linsolve(A,B);
Q4.G = sym2tf(Currents(1)/Vin);

%% Q5
Q5.G = sym2tf(Currents(2)/Vin);