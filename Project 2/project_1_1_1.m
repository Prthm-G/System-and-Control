%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solution to Assignment, Project or Exam
% Populate Q structs with answers (Assig / Exam)
% Populate S structs with answers (Project)
% 
% Once you are done ...
% Run the solution script (this script).
% Run makeMat341.p to generate a MAT file.
%
% Always include a comment above a Matlab function.
% It is automatically turned into a help screen.
%
% Try it. To see this comment:
% In Matlab, enter: help solveAssig
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A blank line terminates the help screen
% This line does not show up in the help screen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize
% Clean up workspace
clear all; clc;
% SN variable must contain Student Number
% This must be right or solution will not be graded
SN    = 18879288;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% student Number
A=9+10;
B=0+10;
C=3+10;
D=0+10;
E=0+10;
F=8+10;
G=1+10;
H=5+10;
s=tf('s');

%---------------Question_1--------------------------%
%figure(1);
% clf;
% p1DSPlot222(SN);
% xlim([0,5]);
% ylim([0,25]);
% xticks(0:1:5);
% yticks(0:1:25);

%units are in volts
Peak=13.6003;
FV=10.995;
Kdc=FV;
%All units are in seconds
Q1_Tr=0.72E-3;
Q1_Tp=1.52E-3;
Q1_Ts=4.6E-3;
Q1_OS=(Peak-FV)/FV;
%all times are in mili seconds
Q1.Tr=Q1_Tr*1E3;
Q1.Tp=Q1_Tp*1E3;
Q1.Ts=Q1_Ts*1E3;
Q1.Pos=Q1_OS*1E2;

%---------------Question_2--------------------------%
%compute the second order approximation that best approximate the raw data
%the approximation should provide a balance between rise an dsettle time
%by using the mean value of the associated natureal frequency
zeta=sqrt((log(Q1_OS))^(2)/((pi^(2))+(log(Q1_OS)^(2))));
beta_rise=sqrt(1-((zeta)^(2)));
wn_rise=(pi-atan(beta_rise/zeta))/(Q1_Tr*beta_rise);
wn_settle=4/(zeta*Q1_Ts);
wn_mean=(wn_rise+wn_settle)/2;
figure(1);
clf;
p1DSPlot222(90300815);
hold on;
%generating the trnasfer function
Q2_rise=(Kdc*(wn_rise)^(2))/((s^(2))+(2*zeta*wn_rise*s)+(wn_rise^(2)));
Q2_settle=(Kdc*(wn_settle)^(2))/((s^(2))+(2*zeta*wn_settle*s)+(wn_settle^(2)));
Q2_mean=(Kdc*(wn_mean)^(2))/((s^(2))+(2*zeta*wn_mean*s)+(wn_mean^(2)));
%generate the data points from the transfer function
[Y1,X1] = step(Q2_rise, 1);
[Y2,X2] = step(Q2_settle, 1);
[Y3,X3] = step(Q2_mean,1);
%plot the graph
%plot(X1*1000,Y1,'r','LineWidth',2);
%plot(X2*1000,Y2,'blue','LineWidth',2);
plot(X3*1000,Y3,'green','LineWidth',2);
%legend("Experimental Curve","Rise time Approximation");
hold off
xlim([0,5]);
ylim([0,25]);
xticks(0:1:5);
yticks(0:1:25);
%legend('Raw curve','Rise','settle','mean');
title("2nd order approximation");
Q2.G=Q2_mean;


%---------------Question_3--------------------------%
G1=A*7/(s+B*800);
G2=C*8/(s+D*700);
G3=10^(5)/(s+E*600);
G4=F*50/(s+G*500);
H1=4/(s+H*5);

G1.u='A1';
G1.y='A2';

G2.u='A3';
G2.y='A4';

G3.u='A5';
G3.y='A6';

G4.u='A6';
G4.y='A7';
H1.u='A8';
H1.y='A9';

sum1 = sumblk('A3 = A2+A6');
sum2 = sumblk('A5 = A1-A9');
sum3 = sumblk('A8 = A4+A7');
sys1=connect(G1,G2,G3,G4,H1,sum1,sum2,sum3,'A1','A4');
Hs=tf(sys1);
%this is in mv/degree 
Q3_TF=((G1*G2*G3*G4*H1)+(G1*G2)+(G2*G3))/((G2*G3*H1)+(G3*G4*H1)+1);
%to get the DC gain we set s=0 notice this is not scaled by QE-3
DC_G=(4.404808136178435e+53/4.274454253365408e+54);
Q3.Kdc=DC_G*1E-3;
%Q3.Kdc=0.1030;
Q3.D=(Q3_TF)/(DC_G);
%Q3.D=Hs/Q3.Kdc;
figure(2);
step(Q3.D)
xlim([0,0.025]);
title("Pure gain of the dynamic sensor");

%---------------Question_4--------------------------%
L6=100*1E-3;%m
Rw=A/2;%ohm
Lw=B*30*1E-6;%H
Km=C*1E-3;%Nm/A
Mm=(D+E)*1E-3;%kg
Bm=(F/30)*1E-6;%Nms
Jr=(G/15)*1E-7;%Kg-m^2
Js=(H/5)*1E-7;%Kg-m^2
Ms=(A/4)*1E-3;%Kg
Bs=(B/3);%Ns/m
Na=(3*1E-2)/(2*pi);%was cm/turn -->m/rad
Jf=3*(C/3)*1E-7;%kg-m^2
Bf=3*D*1E-3;%Nms
Nf=(10*pi/180)/(1E-2);%was degree/cm -->rad/m
Bt=3*E*1E-3;%Nms;
Kt=3*F*30*1E-3;%Nm
Nt=L6;%should be m/rad
Bl=3*(G/5);%Ns/m
Kl=3*(H*15);%N/m


% Vin=1;
% Bsp=Na^(2)*Bs;
% Bfpp=Na^(2)*Nf^(2)*Bf;
% Btpp=Na^(2)*Nf^(2)*Bt;
% Blppp=Na^(2)*Nf^(2)*Nt^(2)*Bl;
% Ktpp=Na^(2)*Nf^(2)*Kt;
% Klppp=Na^(2)*Nf^(2)*Nt^(2)*Kl;
% Mnew=Jr+Js+(Na^(2)*(Mm+Ms))+(Na^(2)*Nf^(2)*Jf);
% 
% %3 by 3 state matrix States are Wr Fkt and Iw
% Q4_A=[((-(Bm+Bsp+Bfpp+Btpp)/Mnew)+((Btpp/Mnew)*(Btpp/(Blppp+Btpp+Klppp/s)))), ((-1/Mnew)+((Btpp/Mnew)*(1/(Blppp+Btpp+(Klppp/s))))),Km/Mnew;
%       Ktpp-((Ktpp*Btpp)/(Blppp+Btpp+(Klppp/s))),(-Ktpp/(Blppp+Btpp+(Klppp/s))), 0;
%      -Km/Lw 0 -Rw/Lw];
% Q4_B=[0;0;1/Lw];
% Q4_C=[0,0,Km;
%       Na/s,0,0;
%       (Na*Nf)/s,0,0;
%       0,0,Km/(3*Na*Nf*Nt);
%       1/s,0,0];
% Q4_D=[0;0;0;0;0];
% 
% %phi= inv(s*I-A);%X=phi*B*u;%Y=(C*Phi*B+D)*u
% %Y=G*u where G=C*Phi*B+D
% Phi=inv(s*eye(3)-Q4_A);
% Q4_Y=(Q4_C*Phi*Q4_B+Q4_D)*Vin;
% Q4_G=Q4_C*Phi*Q4_B+Q4_D;
% Q4.G=Q4_G(1);
% step(Q4.G);
% Q5.G=Q4_G(2);
% step(Q5.G);
% Q6.G=Q4_G(3);
% step(Q6.G);
% Q7.G=Q4_G(4);
% step(Q7.G);
% Q8.G=Q4_G(5)*(180/pi)*Q3.D*Q3.Kdc;
% step(Q8.G);
% 
% % %this Y out put will be torque M
% % sys4=ss(Q4_A,Q4_B,Q4_C,Q4_D);
% % Q4.G=tf(sys4);
Cnew=Jr+Js+(Na^(2)*(Mm+Ms))+(Na^(2)*Nf^(2)*Jf);
Rtt=1/((Bm)+((Na^(2)*Bs))+((Na^(2)*Nf^(2)*Bf)));
Rtpp=1/(Na^(2)*Nf^(2)*Bt);
Ltpp=1/(Na^(2)*Nf^(2)*Kt);
Rippp=1/(Na^(2)*Nf^(2)*Nt^(2)*Bl);
Lippp=1/(Na^(2)*Nf^(2)*Nt^(2)*Kl);
% Zcnew=1/(Cnew*s);
% % Zltpp=s/(Ltpp);
% % Zlippp=s/(Lippp);
% Zltpp=Ltpp*s;
% Zlippp=Lippp*s;
% MX_G=[(1/Zcnew)+(1/Rtt)+(1/Rtpp)+(1/Zltpp),-(1/Rtpp)-(1/Zltpp);-(1/Rtpp)-(1/Zltpp),(1/Rtpp)+(1/Zltpp)+(1/Rippp)+(1/Zlippp)];
% MX_I=[1;0];
% MX_V=MX_G\MX_I;
% Q5_da=Na*MX_V(1)/s;
% step(Q5_da,2000000);
% Q5_qf=Na*Nf*(MX_V(1)/s);
% step(Q5_qf,20000);
% Q5_ft=tf((1/3)*((MX_V(1)-MX_V(2))/(Ltpp*s)));
% step(Q5_ft,2000000);
Q5_A=[((-1/(Cnew*Rtt))+(-1/(Cnew*Rtpp))+((1/(Cnew*Rtpp*Rtpp))/((1/Rippp)+(1/Rtpp)))),((-1/Cnew)+((1/(Cnew*Rtpp))/((1/Rippp)+(1/Rtpp)))),((-1/(Cnew*Rtpp))/((1/Rippp)+(1/Rtpp))),Km/Cnew;
    ((1/(Ltpp))-((1/(Ltpp*Rtpp))/((1/Rippp)+(1/Rtpp)))),(-((1/Ltpp)/((1/Rippp)+(1/Rtpp)))),((1/Ltpp)/((1/Rippp)+(1/Rtpp))),0;
    ((1/(Lippp*Rtpp))/((1/Rippp)+(1/Rtpp))),((1/Lippp)/((1/Rippp)+(1/Rtpp))),(-(1/Lippp)/((1/Rippp)+(1/Rtpp))),0;
    (-Km/Lw),0,0,(-Rw/Lw)];
Q5_B=[0;0;0;1/Lw];
Q5_C=[0,0,0,Km;
      Na,0,0,0;
      Na*Nf,0,0,0
      0,1/(3*Na*Nf*Nt),0,0;];
Q5_D=[0;0;0;0];

sys5=ss(Q5_A,Q5_B,Q5_C,Q5_D)
Q5_TF=tf(sys5)
Q4.G=Q5_TF(1) 
figure(3);
step(Q5_TF(1));
title("Motor Torque Nm Vs input voltage V");

%---------------Question_5--------------------------%
%The Q5_TF (2) is Va to get da we divide by s
Q5.G=(Q5_TF(2)/s);
%Q5.GH=(Q5_TF(2)/s);
figure(4);
step(Q5_TF(2)/s);
title("Actuator displacement m vs Input voltage");

%---------------Question_6--------------------------%
%convert back to degrees
%this is finger angle
Q6.G=(Q5_TF(3)/s)*(180/pi);
figure(5);
step((Q5_TF(3)/s)*(180/pi));
title("Finger angle degree VS input voltage");

%---------------Question_7--------------------------%
%tine force
Q7_A=[((-1/(Cnew*Rtt))+(-1/(Cnew*Rtpp))+((1/(Cnew*Rtpp*Rtpp))/((1/Rippp)+(1/Rtpp)))),((-1/Cnew)+((1/(Cnew*Rtpp))/((1/Rippp)+(1/Rtpp)))),((-1/(Cnew*Rtpp))/((1/Rippp)+(1/Rtpp))),Km/Cnew;
    ((1/(Ltpp))-((1/(Ltpp*Rtpp))/((1/Rippp)+(1/Rtpp)))),(-((1/Ltpp)/((1/Rippp)+(1/Rtpp)))),((1/Ltpp)/((1/Rippp)+(1/Rtpp))),0;
    ((1/(Lippp*Rtpp))/((1/Rippp)+(1/Rtpp))),((1/Lippp)/((1/Rippp)+(1/Rtpp))),(-(1/Lippp)/((1/Rippp)+(1/Rtpp))),0;
    (-Km/Lw),0,0,(-Rw/Lw)];
Q7_B=[0;0;0;1/Lw];
Q7_C=[1,0,0,0;
      (1/Rtpp)/((1/Rippp)+(1/Rtpp)),1/((1/Rippp)+(1/Rtpp)),-1/((1/Rippp)+(1/Rtpp)),0;];
Q7_D=[0;0;];
sys7=ss(Q7_A,Q7_B,Q7_C,Q7_D);
Q7_TF=tf(sys7);

Q7_V1=Q7_TF(1);
Q7_V2=Q7_TF(2);
Q7_Rizz=((Q7_V1-Q7_V2)/(Rtpp))/(3*Na*Nf*Nt);
Q7.G=Q7_Rizz+Q5_TF(4);
figure(6);
step(Q7.G);
title("Tine Force Input force N vs Tine");

%---------------Question_8--------------------------%
%wants degree as input and voltage as output
Q8.GH=((Q5_TF(2)*(1/(Na)))/s)*(180/pi)*Q3.D*Q3.Kdc*Q2.G;
figure(7);
step(Q8.GH);
title("Sensor signal VS AMP voltage");
