function dynamics = getDynamicsPID(Kp, Ki, Kd, CF, Nhat)
    s = tf('s');
    dynamics = Kp + Ki * (1/s) + Kd * ((2*CF*s) / (Nhat*s + 2*CF));
end