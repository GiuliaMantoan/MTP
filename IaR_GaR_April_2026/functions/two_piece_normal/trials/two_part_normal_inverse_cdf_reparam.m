function x = two_part_normal_inverse_cdf_reparam(F, mode, lambda1, lambda2)
    % Inverse CDF of a two-piece normal, param'd by (mode, lambda1, lambda2).
    % F may be a vector of probabilities in [0,1].

    % 1) turn logs back into positive sigmas
    sigma1 = exp(lambda1);
    sigma2 = exp(lambda2);

    % 2) compute cutoff = P(X <= mode)
    cutoff = sigma1 / (sigma1 + sigma2);

    % 3) sanity check
    if any(F < 0) || any(F > 1)
        error('F must be in [0,1].');
    end

    % 4) allocate output
    x = nan(size(F));

    % 5) left branch  (F <= cutoff)
    iL = (F <= cutoff);
    if any(iL)
        % normalise F into [0,1] on the left piece
        uL = F(iL) / cutoff;
        x(iL) = mode + sigma1 * norminv(uL);
    end

    % 6) right branch (F > cutoff)
    iR = (F > cutoff);
    if any(iR)
        % subtract off the left mass, then renormalise
        uR = (F(iR) - cutoff) / (1 - cutoff);
        x(iR) = mode + sigma2 * norminv(uR);
    end
end