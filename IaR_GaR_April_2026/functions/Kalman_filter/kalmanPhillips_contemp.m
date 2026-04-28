function [nLL,innov,x_filt,P_filt,x_smooth,P_smooth] = ...
    kalmanPhillips_contemp(theta,y,drv,ratioVar)
% Measurement gap uses  u_{t} - u^*_t  & oil terms Dz_t, Dz_{t-1}, Dz_{t-2}
% θ = [α1 α2 β γ0 γ1 γ2 logVar_eps x0]

a1   = theta(1);
a2   = theta(2);
beta = theta(3);
g0   = theta(4);
g1 = theta(5);
g2 = theta(6);
var_eps = exp(theta(7));
var_eta = ratioVar*var_eps;
x_prev  = theta(8);  % initial guess
P_prev = 5;

dpi_l1 = drv(:,1);
dpi_l2 = drv(:,2);
dz_0   = drv(:,3);
dz_l1 = drv(:,4);
dz_l2 = drv(:,5);
u_0    = drv(:,6);

T  = numel(y);
x_filt = zeros(T,1);
P_filt = zeros(T,1);
x_pred = zeros(T,1);
P_pred = zeros(T,1);
innov  = zeros(T,1);
S_t    = zeros(T,1);
nLL = 0;

for t = 1:T
    % prediction
    x_pred(t) = x_prev;
    P_pred(t) = P_prev + var_eta;

    % deterministic part of measurement (note −β*u_t)
    c_t = a1*dpi_l1(t) + a2*dpi_l2(t) - beta*u_0(t) + ...
        g0*dz_0(t) + g1*dz_l1(t) + g2*dz_l2(t);

    % innovation
    S = beta^2 * P_pred(t) + var_eps;
    K = (P_pred(t)*beta) / S;
    v = y(t) - c_t - beta*x_pred(t);

    % update
    x_filt(t) = x_pred(t) + K*v;
    P_filt(t) = (1-K*beta)*P_pred(t);

    nLL = nLL + 0.5*(log(2*pi)+log(S)+(v^2)/S);

    x_prev = x_filt(t);  P_prev = P_filt(t);
    innov(t)=v;  S_t(t)=S;
end

% Rauch-Tung-Striebel smoother
x_smooth = x_filt;   P_smooth = P_filt;
for t=T-1:-1:1
    C = P_filt(t)/P_pred(t+1);
    x_smooth(t) = x_filt(t) + C*(x_smooth(t+1)-x_pred(t+1));
    P_smooth(t) = P_filt(t) + C^2*(P_smooth(t+1)-P_pred(t+1));
end
end