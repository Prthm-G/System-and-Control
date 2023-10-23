function Kd = getKdPID(z1_neg, z2_neg, p1_neg, Nhat, Kp)
    Kd = (1/(z1_neg*z2_neg)) + (Kp/(p1_neg/Nhat));
end