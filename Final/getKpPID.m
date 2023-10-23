function Kp = getKpPID(z1_neg, z2_neg, p1_neg, Nhat)
    Kp = 1/(p1_neg/Nhat) - ((z1_neg+z2_neg)/(z1_neg*z2_neg));
end