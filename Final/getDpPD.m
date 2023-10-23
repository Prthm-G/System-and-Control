function dynamics = getDpPD(p1_neg, Nhat)
    s = tf('s');
    dynamics = (-p1_neg/Nhat) * (1/(s-(p1_neg/Nhat)));
end