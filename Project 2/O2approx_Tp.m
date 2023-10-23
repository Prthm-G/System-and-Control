%O2approx_Tp.m
%Generates a second order approzimation given Kdc, os, and Tp
%os in decimal (ex: 20% -> 0.2)
%syntax:
%   [Z,Wn,G] = O2approx_Tp(Kdc,os,Tp)

function [Z,Wn,xfer] = O2approx_Tp(Kdc,os,Tp)

numCoDamp = (log(os))^2;
denCoDamp = pi^2 + (log(os))^2;
CoDamp = sqrt(numCoDamp/denCoDamp);
Z = CoDamp;

nFreqTp = pi / (sqrt(1 - CoDamp^2) * Tp);
Wn = nFreqTp;

num = [nFreqTp^2];
den = [1 2*CoDamp*nFreqTp nFreqTp^2];
xfer = Kdc * tf(num, den);
end