function [ val ] = qs( q,a,y )
%qs: Vector of Quantile Scores as in Gneiting and Ranjan, JBES, (2011)

%   Feb 2014 Version
%   Author: Michael Smith, University of Melbourne.

%   Inputs a=(a(1),....,a(J)), q=(q(1),..,q(J), y=scalar
%   Output val=(QS(1),...,QS(J)) where
%
%   QS(j)= 2 * [ I(y<q(j)) - a(j)] * (q(j) - y)

J=size(q,1);
if(size(a,1)~=J)
    disp('q and a have to be same size');
    stop;
end;

idx=zeros(J,1);
idx(q > y)=1;
val = 2*(idx - a).*(q-ones(J,1)*y);

end

