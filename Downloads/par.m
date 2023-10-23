% %%%%%%%%%%%%%%%%
%  Function: par
%  Calculates the equivalent impedance of parallel impedances
%
% Author: Brandon Just
% %%%%%%%%%%%%%%%%
function Zeq = par(Zs)
    Zeq = 1/sum(1./Zs);
end