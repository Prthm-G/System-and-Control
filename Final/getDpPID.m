function dynamics = getDpPID(p1_neg, Nhat)
    s = tf('s');
    dynamics = (-p1_neg/Nhat) * (1/(s*(s-(p1_neg/Nhat))));
end