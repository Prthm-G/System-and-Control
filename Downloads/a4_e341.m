% Note: All variables are cleared when this is run
clear all; clc;

%set student number variables
numA = 10 + 1;
numB = 10 + 2;
numC = 10 + 3;
numD = 10 + 4;
numE = 10 + 5;
numF = 10 + 6;
numG = 10 + 7;
numH = 10 + 8;

%Set variables
s = tf('s');

G1 = 1/(s+numA);
G2 = 10/(s+numB);
G3 = 10/(s+numC);
G4 = 10/(s+numD);
G5 = 1/(s+numE);

H1 = 100/(s+numF);
H2 = 1000/(s+numG);

%%%%%%%%%%%
% Q1
%%%%%%%%%%%

if 1
   Q1.G11 = minreal(((G1*G2) / (1+((G1*G2)*(H2*G5)*((G3/G2)+(H1*G4))))))
   Q1.G12 = minreal(((H1*G2*G4 + G3)*(G1)*(G5))/(1+(H1*G2*G4 + G3)*(G1)*(G5)*(H2)));
   
   AG21 = G1*(-1*H2*G5);
   BG21 = AG21/(1-AG21*G3);
   CG21 = BG21*G2*G4;
   DG21 = CG21/(1-CG21*H1);
   Q1.G21 = minreal(-1*DG21);
   
   AG22 = G1*(-1*H2);
   BG22 = G2*G4*H1 + G3;
   CG22 = -1*G4;
   DG22 = G5/(1-(G5*AG22*BG22));
   Q1.G22 = minreal(CG22*DG22);
end

%%%%%%%%%%%
% Q2
%%%%%%%%%%%

if 1
    Q2.G1 = Q1.G21 + Q1.G11;
    Q2.G2 = Q1.G22 + Q1.G12;
end

%%%%%%%%%%%
% Q3
%%%%%%%%%%%

if 1
    %set variables
    Rw = 1 + numA/10;
    Lw = (100 + 10 * numB) * 10^(-6);

    Jr = (numC / 10) * 10^(-6);
    Br=(numD + numE + numF) * 10^(-6);

    km = (10 + numG)/1000;
    Km = (10 + numG)/1000;

    Vin = 12; %input of 1 (for some reson i cant use imput 12/s then divide output by 12/s)

    Jf = numG/30 * 10^(-6);
    Bf = numH * Br;

    tf_elec = minreal(1/(s*Lw + Rw));
    tf_mech = minreal(1/(s*(Jr + Jf) + Br + Bf));

    tf_noFeedback = minreal(tf_elec * km * tf_mech);

    tf_Feedback = minreal(feedback(tf_noFeedback, Km))

    step(tf_Feedback)

    Q3.Ye = tf_elec;
    Q3.Ym = tf_mech;
    Q3.G = tf_Feedback;
end

%%%%%%%%%%%
% MakeMat
%%%%%%%%%%%

if 0
   SN = 12345678
   makeMat341 
end