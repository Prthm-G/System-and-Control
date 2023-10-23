% syntax: [kd, dp] = getDynamicsPD
function [Kd, dynamics] = getDynamicsPD(p_neg, z_neg, Nhat)
    s = tf('s');
    Kd = (z_neg - (p_neg/Nhat)) / (z_neg*(p_neg/Nhat));
    dynamics = ((p_neg/Nhat)/z_neg) * ((s-z_neg)/(s-p_neg));
end