function [ val ] = qwps( x,xobs,J,QW )
%qwps: compute the Quantile Weighted Probabiltiy Score at point xobs, using iterates x from the predictive distribution

%   Feb 2014 Version
%   Author: Michael Smith, University of Melbourne.

%   The formula is from Gneiting and Ranjan (2011; Eqn.17) and is an
%   (accurate) approximation to the quantile weighted continuous ranked
%   probability score in G&R (2011; Eqn.15). 
%  
%       val = (1/J-1) * sum_{j=1}^{J-1}[ v(a(j))*QS(q(j),xobs) ]
%
%   where:
%       a(j) = j/J
%       q(j) = F^{-1}(a(j)) 
%       F is the EDF built from x
%       QS(q,xobs) = quantile score (separate function)
%       v(q) = quantile weight kernel that takes values
%            = 1            if QW=1     (uniform)
%            = q(1-q)       if QW=2     (centre)
%            = (2q-1)^2     if QW=3     (tails)
%            = q^2          if QW=4     (right tail)
%            = (1-q)^2      if QW=5     (left tail)
%            = (2q-1)^4     if QW=6     (heavy tails)


step=1/J; a1=step; a2=1-a1;
a=(a1:step:a2)';   %   a(1),...,a(J)

switch QW       %   compute weighting function values v(1),...,v(J) 
    case 1
        v=ones(J-1,1);
    case 2
        v=a.*(ones(J-1,1)-a);
    case 3
        v=(2.*a-1).^2;
    case 4
        v=a.^2;
    case 5
        v=(1-a).^2;
    case 6
        v=(2.*a-1).^4;
    otherwise
        disp('QW out of range');
end

%% jitter observtions if there are ties to overcome CDF routine constraints
x=sort(x);
x = x(:); 
x = x(~isnan(x)); 

xdiff=diff(x);
eps=mean(xdiff)*0.0001;
k=find(xdiff<=eps);
if(~isempty(k))   % jitter observations by eps        
    jitter=normrnd(0,10*eps,size(k));
    x(k) = x(k) + jitter;
end
%%
%obj=paretotails(xk,0,1,'ecdf');
obj=paretotails(x,0,1,@ecdf);  %fit EDF to the data
q=icdf(obj,a);  %   compute q(1),...,q(J)
qsval=qs(q,a,xobs); %   compute vector of QS(q(1),xobs),...,QS(q(J),xobs)\

val = sum(qsval.*v)/(J-1);

end

