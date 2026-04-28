function [nLL,innov,x_filt,P_filt,x_smooth,P_smooth] = ...
    kalmanPhillips_lag(theta,y,drv,ratioVar)
% Measurement gap uses  u_{t-1} - u^*_{t-1}  (as in original code)
% θ = [α1 α2 β γ1 γ2 logVar_eps x0]

a1   = theta(1);
a2   = theta(2);
beta = theta(3);
g1   = theta(4);
g2   = theta(5);
var_eps = exp(theta(6));
var_eta = ratioVar*var_eps;
x_prev  = theta(7);
P_prev = 5;

dpi_l1 = drv(:,1);
dpi_l2 = drv(:,2);
dz_l1  = drv(:,3);
dz_l2  = drv(:,4);
u_l1   = drv(:,5);

T  = numel(y);
x_filt = zeros(T,1);   
P_filt = zeros(T,1);
x_pred = zeros(T,1);   
P_pred = zeros(T,1);
innov  = zeros(T,1);   
S_t    = zeros(T,1);
nLL = 0;

for t = 1:T
    % prediction (random walk)
    x_pred(t) = x_prev;
    P_pred(t) = P_prev + var_eta;

    % deterministic part of measurement
    c_t = a1*dpi_l1(t) + a2*dpi_l2(t) - beta*u_l1(t) + g1*dz_l1(t) + g2*dz_l2(t);

    % innovation
    S = beta^2 * P_prev + var_eps;                 % note P_prev
    K = (P_prev*beta) / S;
    v = y(t) - c_t - beta*x_prev;                  % uses x_prev ≡ u^*_{t-1}

    % update
    x_filt(t) = x_prev + K*v;
    P_filt(t) = (1-K*beta)*P_prev;

    nLL = nLL + 0.5*(log(2*pi)+log(S)+(v^2)/S);

    x_prev = x_filt(t);
    P_prev = P_filt(t);
    innov(t)=v;  S_t(t)=S;
end

% Rauch-Tung-Striebel smoother
x_smooth = x_filt;
P_smooth = P_filt;

for t=T-1:-1:1
    C = P_filt(t)/P_pred(t+1);
    x_smooth(t) = x_filt(t) + C*(x_smooth(t+1)-x_pred(t+1));
    P_smooth(t) = P_filt(t) + C^2*(P_smooth(t+1)-P_pred(t+1));
end
end