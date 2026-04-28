function [mean_value, std_dev, skewness, kurtosis, delta, mu,omega] = moments_skewt_distr_sm_updated(loc, scale, shape, df)
    
    % Step 2: Calculate delta
    delta = shape / sqrt(1 + (shape^2));
    
    % Step 3: Calculate mu
    mu = delta * sqrt(df/pi) * (gamma((0.5 * (df - 1))) / gamma((0.5 * df)));
    
    % Step 4: Set omega equal to scale
    omega = scale; 
    
    % Step 5: Calculate mean
    mean_value = loc + omega * mu;
    
    % Step 6: Calculate the second moment
    second_moment = (loc^2) + 2 * loc * omega * mu + ((omega^2 * df) / (df - 2));
    
    % Step 7: Calculate variance
    variance = second_moment - (mean_value^2);
    
    % Step 7.1: Calculate sd
    std_dev = sqrt(variance);

    % Step 8: Calculate the third moment
    third_moment = (loc^3) + 3 * (loc^2) * omega * mu + 3 * loc * ((omega^2 * df) / (df - 2)) + ...
                   omega^3 * mu * (3 - (delta^2)) * df / (df - 3);

    % Step 9: Calculate the third central moment
    third_central_moment = third_moment - 3 * mean_value * second_moment + 2 * (mean_value^3);
    
    % Step 10: Calculate skewness
    skewness = third_central_moment / (std_dev^3);

    % Step 11: calculate the fourth moment 
    fourth_moment = (loc^4) + 4* (loc^3) * omega * mu + 6 * (loc^2) *  ((omega^2 * df) / (df - 2)) + 4 * loc *  (omega^3) * mu * (3 - (delta^2)) * df / (df - 3) + (omega ^4) * 3 * (df^2) / ((df - 2) * (df -4));

    % step 12: calculate the third central moment 
    fourth_central_moment = fourth_moment - 4 * third_moment * mean_value + 6 * second_moment * (mean_value^2) - 3 * (mean_value^4);
                    
    % step 13: calculate kurtosis
    kurtosis = (fourth_central_moment / (std_dev^4));

end