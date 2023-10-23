function [Kp, dynamics] = getDynamicsPI(z_neg)
    s = tf('s');
    Kp = -1/z_neg;
    dynamics = 1/s + Kp;
end