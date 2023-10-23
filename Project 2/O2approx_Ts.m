%O2approx_Ts.m
%Generates a second order approzimation given Kdc, os, and Ts
%os in decimal (ex: 20% -> 0.2)
%syntax:
%   [Z,Wn,G] = O2approx_Ts(Kdc,os,Ts)

function [Z,Wn,xfer] = O2approx_Ts(Kdc,os,Ts)

numCoDamp = (log(os))^2;
denCoDamp = pi^2 + (log(os))^2;
CoDamp = sqrt(numCoDamp/denCoDamp);
Z = CoDamp;

nFreqTs = (4 / (CoDamp * Ts));
Wn = nFreqTs;

num = [nFreqTs^2];
den = [1 2*CoDamp*nFreqTs nFreqTs^2];
xfer = Kdc * tf(num, den);
end