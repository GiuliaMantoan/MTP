function [beta, beta_ols, se, se_boot] = quant_reg(X, y, tau, stderror)

    % X: Matrix of predictors (n x p)
    % y: Response vector (n x 1)
    % tau: Quantile level (e.g., 0.5 for median regression)

    % Ensure the inputs are column vectors/matrices
    y = y(:);
    [n, p] = size(X);
    
    % Add an intercept term to X
    X = [ones(n, 1), X];
    
    % Initial estimates using linear regression
    beta_ols = (X' * X) \ (X' * y);
    
    % Define the quantile loss function
    quantile_loss = @(beta) sum((tau - (y - X * beta < 0)) .* (y - X * beta));
    
    % Use fminsearch to minimize the quantile loss function
    options = optimset('Display', 'off', 'TolX', 1e-8, 'TolFun', 1e-8);
    beta = fminsearch(quantile_loss, beta_ols, options);
    
    % Initialize outputs
    se = [];
    se_boot = [];
        
    % If standard errors are requested
    if stderror == 1
            
        % fitted values 
        yfit= X * beta;
    
        % Calculate residuals
        residuals = y - yfit;
    
        % Calculate the asymptotic covariance matrix
        rho_tau = tau - (residuals < 0);
        X_tilde = X .* rho_tau;
        H = (X' * X) / n;
        Omega = (X_tilde' * X_tilde) / n;
    
    
        % Calculate standard errors
        V = H \ Omega / H / n;
        se = sqrt(diag(V));
    
    
        % bootstrapped std error
        rho=@(r)sum(abs(r.*(tau-(r<0)))); % quantile loss function 
        %pboot=bootstrp(200,@(bootr)fminsearch(@(beta)rho(yfit+bootr-X * beta),beta)', residuals);
        pboot = bootstrp(200, @(bootr)fminsearch(@(beta)rho(yfit + bootr - X * beta), se, options)', residuals);
    
        % calculate bootstrapped se 
        se_boot=std(pboot);

    end 
    
end
