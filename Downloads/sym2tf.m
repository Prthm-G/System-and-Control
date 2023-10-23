% %%%%%%%%%%%%%%%%
%  Function: sym2tf
%  Converts a symbolic expression to a transfer function
%  Author: Brandon Just
% %%%%%%%%%%%%%%%%
function TF = sym2tf(sym)
    [num, den] = numden(sym);
    tfnum = sym2poly(num);
    tfden = sym2poly(den);
    TF = minreal(tf(tfnum, tfden));
end