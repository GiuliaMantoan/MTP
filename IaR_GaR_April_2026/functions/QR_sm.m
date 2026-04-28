function Y_f=QR_sm(Q_fit, quantilelevels)

% INPUT
% X: explanatory variables
% y: dep variables
% quantilelevels = quantiles 

%[T,p]=size(X);
k=length(quantilelevels); % # of quantiles
T = size(Q_fit,1);

%% STEP1 : Run the quantile regressions

%beta=zeros(p,k); % pre-allocate space for coeff

%for r=1:k

%   beta(:,r)=rq(X, y, quantilelevels(r)); %loop over each quantile 

%end  

%%  STEP 2: In-Sample density fit

%Q_fit=X*beta; % get the fitted quantiles of the distr 
Q_fit=sort(Q_fit,2); % sort the quantiles 
N=20000;% total sample size 
n=N/(k+1);

Y_f=zeros(T,N);

for j=2:k
    ind=(j-1)*n+1:j*n; % get equal space within quantile 2 up to quantile end-1
    Y_f(:,ind)=Q_fit(:,j-1)+sparse(1:T,1:T,Q_fit(:,j)-Q_fit(:,j-1))*rand(T,n); % q_{i-1} + (q_{i} -q_{i-1} * random uniform that covers the equal space interval
end
sig_l=(Q_fit(:,2)-Q_fit(:,1))./(norminv(quantilelevels(2))-norminv(quantilelevels(1)));
mu_l=Q_fit(:,1)-sig_l*norminv(quantilelevels(1));
Y_f(:,1:n)=mu_l+sparse(1:T,1:T,sig_l)*norminv(quantilelevels(1)*rand(T,n));
sig_u=(Q_fit(:,k)-Q_fit(:,k-1))./(norminv(quantilelevels(k))-norminv(quantilelevels(k-1)));
mu_u=Q_fit(:,k)-sig_u*norminv(quantilelevels(k));
Y_f(:,k*n+1:end)=mu_u+sparse(1:T,1:T,sig_u)*norminv(quantilelevels(k)+(1-quantilelevels(k))*rand(T,n));

N = size(Y_f, 2);
