function gam = stdtogam(mean_val, mod_val, std_val)

    % Compute the intermediate variable c (this is what in wikipedia is
    % called beta
    c = (((mean_val - mod_val)^2) * pi) / (2 * (std_val^2));
    
    % Conditional computation of gam
    if (mean_val - mod_val) < 0
        gam = -sqrt(1 - (((sqrt(1 + 2 * c) - 1) / c)^2));
    else
        gam = sqrt(1 - (((sqrt(1 + 2 * c) - 1) / c)^2));
    end
end