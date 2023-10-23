function [dynamics, zeros, max_PM] = getComplexZeros(crossoverMagnitude, Kref, poleDynamics, openLoop)
    s = tf('s');
    direction = 1;
    vertical = 1;
    horizontal = 1;
    max_PM = 0;
    prev_PM = 0;
    z1_test = intmin;
    z2_test = intmin;
    z1_res = intmin;
    z2_res = intmin;
    D_res = s;
    resolution = 0.001;
    sigma = -crossoverMagnitude;
    omega = 0;
    
    for iter = 1:1000
        if direction == 1
            if vertical == 1
                while 1
                    z1_test = sigma + omega*1i;
                    z2_test = sigma - omega*1i;
                    test_Dz = getDzPID(z1_test, z2_test);
                    [~, test_PM] = margin(Kref * poleDynamics * test_Dz * openLoop);
                    if (test_PM) < (prev_PM)
                        omega = omega - resolution;
                        if direction == 1
                            direction = 0;
                        else
                            direction = 1;
                        end
    
                        if vertical == 1
                            vertical = 0;
                        else
                            vertical = 1;
                        end
                        break
                    else
                        if (test_PM) > (max_PM)
                            max_PM = test_PM;
                            z1_res = z1_test;
                            z2_res = z2_test;
                            D_res = test_Dz * poleDynamics;
                            prev_PM = test_PM;
                        end
                    end
                    omega = omega + resolution;
                end
            else
                while 1
                    if omega - resolution < 0
                        break
                    end
                    z1_test = sigma + omega*1i;
                    z2_test = sigma - omega*1i;
                    test_Dz = getDzPID(z1_test, z2_test);
                    [~, test_PM] = margin(Kref * poleDynamics * test_Dz * openLoop);
                    if (test_PM) < (prev_PM)
                        omega = omega + resolution;
                        if direction == 1
                            direction = 0;
                        else
                            direction = 1;
                        end
    
                        if vertical == 1
                            vertical = 0;
                        else
                            vertical = 1;
                        end
                        break
                    else
                        if (test_PM) > (max_PM)
                            max_PM = test_PM;
                            z1_res = z1_test;
                            z2_res = z2_test;
                            D_res = test_Dz * poleDynamics;
                            prev_PM = test_PM;
                        end
                    end
                    omega = omega - resolution;
                end
            end
        else
            if horizontal == 1
                while 1
                    if sigma + resolution >= 0
                        break;
                    end
                    z1_test = sigma + omega*1i;
                    z2_test = sigma - omega*1i;
                    test_Dz = getDzPID(z1_test, z2_test);
                    [~, test_PM] = margin(Kref * poleDynamics * test_Dz * openLoop);
                    if (test_PM) < (prev_PM)
                        sigma = sigma - resolution;
                        if direction == 1
                            direction = 0;
                        else
                            direction = 1;
                        end
    
                        if horizontal == 1
                            horizontal = 0;
                        else
                            horizontal = 1;
                        end
                        break
                    else
                        if (test_PM) > (max_PM)
                            max_PM = test_PM;
                            z1_res = z1_test;
                            z2_res = z2_test;
                            D_res = test_Dz * poleDynamics;
                            prev_PM = test_PM;
                        end
                    end
                    sigma = sigma + resolution;
                end
            else
                while 1
                    z1_test = sigma + omega*1i;
                    z2_test = sigma - omega*1i;
                    test_Dz = getDzPID(z1_test, z2_test);
                    [~, test_PM] = margin(Kref * poleDynamics * test_Dz * openLoop);
                    if (test_PM) < (prev_PM)
                        sigma = sigma + resolution;
                        if direction == 1
                            direction = 0;
                        else
                            direction = 1;
                        end
    
                        if horizontal == 1
                            horizontal = 0;
                        else
                            horizontal = 1;
                        end
                        break
                    else
                        if (test_PM) > (max_PM)
                            max_PM = test_PM;
                            z1_res = z1_test;
                            z2_res = z2_test;
                            D_res = test_Dz * poleDynamics;
                            prev_PM = test_PM;
                        end
                    end
                    sigma = sigma - resolution;
                end
            end
        end

        dynamics = D_res;
        zeros = [z1_res, z2_res];
    end