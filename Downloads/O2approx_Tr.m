%O2approx_Tr.m
%Generates a second order approzimation given Kdc, os, and Tr
%os in decimal (ex: 20% -> 0.2)
%syntax:
%   [Z,Wn,G] = O2approx_Tr(Kdc,os,Tr)

function [Z,Wn,xfer] = O2approx_Tr(Kdc,os,Tr)

numCoDamp = (log(os))^2;
denCoDamp = pi^2 + (log(os))^2;
CoDamp = sqrt(numCoDamp/denCoDamp);
Z = CoDamp;

nFreqTr = (1 / (sqrt(1 - CoDamp^2) * Tr)) * (pi - atan(sqrt(1 - CoDamp^2) / CoDamp));
Wn = nFreqTr;

num = [nFreqTr^2];
den = [1 2*CoDamp*nFreqTr nFreqTr^2];
xfer = Kdc * tf(num, den);
end