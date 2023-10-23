function K = getPhaseMarginK(Kref, fullDynamics, openLoop, targetPhaseMargin)
    K = Kref;
    while 1
        [~, test_PM] = margin(K * fullDynamics * openLoop);
        if abs(test_PM) <= targetPhaseMargin
            break
        end
        K = K + 0.001;
    end
end

