% Note: All variables are cleared when this is run
clear all; clc;

%set student number variables
numSTU = 90300815;
Ndig = dec2base(numSTU,10) - '0';
SN=90300815;

numA = 10 + Ndig(1);
numB = 10 + Ndig(2);
numC = 10 + Ndig(3);
numD = 10 + Ndig(4);
numE = 10 + Ndig(5);
numF = 10 + Ndig(6);
numG = 10 + Ndig(7);
numH = 10 + Ndig(8);

%Set variables
s = tf('s');

%%%%%%%%%%%
% Q1
%%%%%%%%%%%
    
    p1DSPlot222(90300815,1);
    %%stepinfo(p1DSPlot222(90300815))
    Q1.FV = 10.995;
    Q1.Kdc = 10.995;
    Q1.Pv = 13.6003;
    Q1.Pos = (Q1.Pv-Q1.Kdc)/Q1.Kdc * 100;
    Q1.Tp = 1.52;
    Q1.Tr = 0.72;
    Q1.Ts = 4.6;
    
    yline(Q1.FV*0.98,'b');
    yline(Q1.FV*1.02,'b');
    yline(Q1.FV,'r');
    xline(Q1.Tp,'r');
    xline(Q1.Tr,'r');
    xline(Q1.Ts,'r');

    Q1.os = (16.2-Q1.Kdc)/Q1.Kdc;
    Q1.TpS = 1.52/1000;
    Q1.TrS = 0.72/1000;
    Q1.TsS = 4.75/1000;
    

%%%%%%%%%%%
% Q2
%%%%%%%%%%%

    
    [Z_Tr,Wn_Tr,xfer_Tr] = O2approx_Tr(Q1.Kdc,Q1.os,Q1.TrS);
    [Z_Ts,Wn_Ts,xfer_Ts] = O2approx_Ts(Q1.Kdc,Q1.os,Q1.TsS);
    
    Z = (Z_Tr+Z_Ts)/2;
    Wn = (Wn_Tr+Wn_Ts)/2;
    xfer_num = Q1.Kdc * Wn^2;
    xfer_den = s^2 + s*Z*Wn*2 + Wn^2;
    xfer = xfer_num/xfer_den;
    %xfer_real = 8.828e07 / (s^2 + 2104*s + 6.791e06)
    
    figure(2);
    hold on
    step(xfer_Ts,'g')
    step(xfer_Tr,'b')
    step(xfer,'y')
    %step(xfer_real,'r')
    hold off
    
    Q2.G = xfer;
   

%%%%%%%%%%%
% Q3
%%%%%%%%%%%

    G1 = numA*7/(s + numB*800); 
    G2 = numC*8/(s + numD*700);
    G3 = 10^5/(s + numE*600); 
    G4 = numF*50/(s + numG*500);
    H1 = 4/(s + numG*5);
    
    Q3_A = G3 / (1+G3*G4*H1);
    Q3_B = 1+Q3_A/G1;
    Q3_C = 1/(1+Q3_A*G2*H1);
    Q3.Hs = minreal(Q3_B*Q3_C*G1*G2);
    Q3.Hs2 = minreal((G2*G3+G1*G2+G1*G2*G3*G4*H1)/(1+G3*G4*H1+G2*G3*H1));%mV
    Q3.Kdc = margin(Q3.Hs2)
    
    Q3.D=Q3.Hs2/Q3.Kdc;%mV
    
    figure(3)
    step(Q3.Hs);
    figure(31)
    step(Q3.Hs2);
    
    

%%%%%%%%%%%
% Q4
%%%%%%%%%%%
   
    % Motor & Lead-Screw Variables
    Rw = numA / 2;      % (Ω) CHECK
    Lw = numB * 30;     % (μH) CHECK
    Km = numC;          % (mNm/A) CHECK
    Mm = numD + numE;   % (g) CHECK
    Jr = numG / 15;     % (g‐cm2) CHECK
    Bm = numF / 30;     % (μNms) CHECK
    Js = numH / 5;      % (g‐cm2) CHECK
    Ms = numA / 4;      % (g) CHECK
    Bs = numB / 3;      % (Ns/m) CHECK
    Na = 3;
    % (cm/turn) CHECK
    
    % Mechanism & Controller Variables
    Jf = numC / 3;      % (g‐cm2) CHECK
    Bf = numD / 50;     % (Nms) CHECK
    Nf = 10;            % (deg/cm) CHECK
    Bt = numE;          % (mNms) CHECK
    Kt = numF * 30;     % (mNm) CHECK
    L6 = 100;           % (mm) CHECK
    Bl = numG/5;
    Kl = numH*15;
    %CF = 200;           % (Hz) CHECK
    
    % Motor & Lead-Screw SI Variables
    Rw = Rw* 1;
    Lw = Lw * 10^(-6);
    Km = Km * 10^(-3);
    Mm = Mm * 10^(-3);
    Jr = Jr * 10^(-7);
    Bm = Bm * 10^(-6);
    Js = Js * 10^(-7);
    Ms = Ms * 10^(-3);
    Bs = Bs * 1;
    Na = Na * 10^(-2)/(2*pi);
   
   
    
    % Mechanism & Controller SI Variables
    Jf = 3 * Jf * 10^(-7);
    Bf = 3 * Bf * 1;
    Nf = Nf * 10^(2) * (pi/180);
    Bt = 3 * Bt * 10^(-3);
    Kt = 3 * Kt * 10^(-3);
    L6 = L6 * 10^(-3);
    Nt = L6*(pi/180);
    %CF = CF * 1;
    
    % Gravity
    g = 9.81;
    %Nt = (Nf * Ns * P) / (2 * pi * R);
    Bfr = (1*Na^2)*(1*Nf^2)*Bf;
    Btr = (1*Na^2)*(1*Nf^2)*Bt;
    Bsr = (1*Na^2)*Bs;
    Ktr = (1*Na^2)*(1*Nf^2)*Kt;
    Klr = (1*Na^2)*(1*Nf^2)*(1*Nt^2)*Kt;
    Blr = (1*Na^2)*(1*Nf^2)*(1*Nt^2)*Bl;
    Mr = Jr + Js + 1*Na^2*(Mm+Ms)+(1*Na^2)*(1*Nf^2)*Jf;

   
    %sQr = 0 Qr + 1 Wr + 0 Iw + 0 Vin;
    %sWr = (-Ktr* 1/Mr) Qr + ((Br + Bsr + Bfr + Btr) * 1/Mr) Wr + (Km * 1/Mr) Iw + 0 Vin;
    %sIw = 0 Qr + (-Km/Lw) Wr + (-Rw/Lw) Iw + (1/Lw) Vin;
    

    Q4.A = [
        0 1 0;
        (-Ktr * 1/Mr) (-(Bl + Bsr + Bfr + Btr+ Blr)* 1/Mr) (Km * 1/Mr);
        0 (-Km/Lw) (-Rw/Lw);
        ];
    
    Q4.B = [
        0;
        0;
        (1/Lw);
        ];
    
    Q4.C = [
        1 0 0; 
        (Na*Nf*Kt/0.1 / 3) 0 0;
        ];
    
    Q4.D = [
        0;
        0;
        ];
     
    Q4.G = tf(ss(Q4.A, Q4.B, Q4.C, Q4.D));
    
    Q4.Gj = Q4.G(1) * 180/pi;
    Q4.Gt = minreal(Q4.G(2) * 1000 / g);

    %Q4.Gt_real = 1.265e09 / (s^3 + 2.491e04 * s^2 + 1.159e08 * s + 1.642e08)

    figure(4)
    step(Q4.Gj)
    figure(41)
    step(Q4.Gt)


%%%%%%%%%%%
% Q5
%%%%%%%%%%%

    % Q5
    
    % Find Du, Kfb = bottom feedback of condensed block diagram, GH
    % Find open loop transfer function
    % Kfb = margin(GH)
    
    % since we have the closed loop transfer function 
    CF = 200; % Hz
    
    Q5.Du = CF/(1*s+CF); % If we know our 10x most dominant pole and CF, how do we calculate the filter?
    Q5.Kfb = 1000/(dcgain(Q3.Hs))
    
    % a CL transfer func = G/(1+G)
    % we have the transfer func Hs
    
    Q5.G = (Q2.G * Q4.Gj * Q3.Hs/1000 * Q5.Kfb * Q5.Du);


    %Q5.GH_real = (3.104e27*s^8 + 2.632e32*s^7 + 9.515e36*s^6 + 1.901e41*s^5 + 2.268e45*s^4 + 1.618e49*s^3 + 6.406e52*s^2 + 1.104e56*s + 9.921e57)/(s^16 + 1.335e05*s^15 + 8.075e09*s^14 + 2.926e14*s^13 + 7.071e18*s^12 + 1.202e23*s^11 + 1.474e27*s^10 + 1.319e31*s^9 + 8.599e34*s^8 + 4.047e38*s^7 + 1.351e42*s^6+ 3.127e45*s^5 + 4.812e48*s^4 + 4.263e51*s^3 + 9.834e53*s^2 + 6.048e55*s + 8.373e55)

    figure(5)
    step(Q5.G)

    
    % Find feedback gain
    % 2<kfb<kultimate
    % ku = margin(GH)


%%%%%%%%%%%
% Q6
%%%%%%%%%%%

    % Q6 
    % Find Xj
    Q6.G = 1*Q2.G*Q4.Gj;
    Q6.H6 = Q3.Hs * 10^(-3) *Q5.Du*Q5.Kfb;
    Q7.G = Q6.G/(1+Q6.H6*Q6.G);
    figure(6)
    step(Q6.G)

%%%%%%%%%%%
% Q7
%%%%%%%%%%%

    % Q7
    % Find Xt
    Q8.GH = Q7.G / Q4.Gj * Q4.Gt;

    %Q7.Xt_real = (1.116e17*s^13 + 1.467e22*s^12 + 8.698e26*s^11 + 3.073e31*s^10 + 7.188e35*s^9 + 1.169e40*s^8 + 1.35e44*s^7 + 1.108e48*s^6 + 6.348e51*s^5 + 2.428e55*s^4 + 5.659e58*s^3 + 6.494e61*s^2 + 1.577e64*s + 9.715e65)/(s^18 + 1.584e05*s^17 + 1.152e10*s^16 + 5.092e14*s^15 + 1.529e19*s^14 + 3.302e23*s^13 + 5.286e27*s^12 + 6.382e31*s^11 + 5.853e35*s^10 + 4.075e39*s^9 + 2.14e43*s^8 + 8.368e46*s^7 + 2.393e50*s^6 + 4.866e53*s^5 + 6.655e56*s^4 + 5.221e59*s^3 + 1.256e62*s^2 + 2.005e64*s + 1.159e66)
  
    figure(7)
    step(Q7.G)

    % Go to the Simulink model and look through what you need to multiply to
    % get to the end position you want 

%%%%%%%%%%%
% makeMat
%%%%%%%%%%%
%%makeMat341



