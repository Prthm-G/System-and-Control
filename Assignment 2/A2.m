%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solution to Assignment 1 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize

% Clean up workspace
clear all; clc;
SN    = 90300815;
Kdc   = 25.7541;

%   QUESTION 1  % 

a2DSPlot(SN,1);

Q1.Tr  = 138;
Q1.Tr1 = 95;  
Q1.Tau = 200;
Q1.Tp   = 370;
Q1.Ts  = 800;
Q1.Pos = 37.258;

%   QUESTION 2  % 
Q2.zeta = 0.299809473154;

%   QUESTION 3  %
Q3.wn  = 14.2443034128;
Q3.wn1 = 59.4428744;
num3 = Kdc*202.900179;
den3 = [1 8.54115420327  202.900179];
Q3.G = tf(num3,den3);
[ys, ts] = step(Q3.G);

figure(6);
clf;
plot(ts, ys, 'b-', 'LineWidth', 3);
grid on;
title('Step Response');                                  
xlabel('Time (sec)');
ylabel('OP Voltage (V)');


%   QUESTION 4  %

Q4.wn = 8.9002085;
Q4.wn1 = 25.6117543539;
num4 = Kdc*79.21371266;
den4 = [1 5.33673368 79.213726644];
Q4.G = tf(num4,den4);

%   QUESTION 5  %
Q5.wn = 16.677258;
num5 = Kdc*278.130941671;
den5 = [1 9.9999999999 278.130941671];
Q5.G = tf(num5,den5);



%   QUESTION 6  %
a6 = (Q3.wn+ Q4.wn)/2;
Q6.wn = a6;
num6  = Kdc*a6^2;
den6  = [1 2*Q2.zeta*a6 a6^2];
Q6.G  = tf(num6,den6);

[ys, ts] = step(Q6.G);

figure(2);
clf;
plot(ts, ys, 'b-', 'LineWidth', 3);
grid on;
title('Step Response');
xlabel('Time (sec)');
ylabel('OP Voltage (V)');

% set time vector
t = 0:1e-3:1;
% get step response vector
stepResponse = step(Q6.G, t);

% plot on figure 1
a2DSPlot(SN, 3)
hold on
plot(stepResponse, 'LineWidth', 2)

%   QUESTION 7  %
a7 = (Q3.wn + Q5.wn)/2;
Q7.wn = a7;
num7  = 25.7541*a7^2;
den7  = [1 2*Q2.zeta*a7 a7^2];
Q7.G  = tf(num7,den7);

[ys, ts] = step(Q7.G);

figure(4);
clf;
plot(ts, ys, 'b-', 'LineWidth', 3);
grid on;
title('Step Response');                                  
xlabel('Time (sec)');
ylabel('OP Voltage (V)');
t = 0:1e-3:1;
% get step response vector
stepResponse = step(Q7.G, t);

% plot on figure 1
a2DSPlot(SN, 5)
hold on
plot(stepResponse, 'LineWidth', 2)

% Question 8 %

peak_time = 0.226453;
DC_gain   = Kdc ;
Os        = 18.586166828;
Q8.zeta   = 0.472168330703;
Q8.wn     = 15.7378476497;
num8      = Kdc*(Q8.wn)^2;
den8      = [1 2*Q8.zeta*Q8.wn (Q8.wn)^2];
Q8.G      = tf(num8,den8);

[ys, ts] = step(Q8.G);

figure(8);
clf;
plot(ts, ys, 'b-', 'LineWidth', 3);
grid on;
title('Step Response');                                  
xlabel('Time (sec)');
ylabel('OP Voltage (V)');

% Question 9 %

Q9.wn = 5;
num9  = Kdc*(Q9.wn)^2;
den9  = [1 2*Q9.wn (Q9.wn)^2];
Q9.G  = tf(num9,den9);
[ys, ts] = step(Q9.G);

figure(9);
clf;
plot(ts, ys, 'b-', 'LineWidth', 3);
grid on;
title('Step Response');                                  
xlabel('Time (sec)');
ylabel('OP Voltage (V)');








if 0
end