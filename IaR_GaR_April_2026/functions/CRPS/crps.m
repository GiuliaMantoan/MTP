function [ val ] = crps( x,xobs,method );
%crps: Compute the CRPS of Gneiting and Raftery (2007) at point xobs, using iterates x from the predictive distribution
% 	Feb 2014 Version
% 	Author: Michael Smith, Univeristy of Melbourne.
%	Uses numerical integration, as in Appendix B of Smith and Vahey (2015) 'Asymmetric Forecast Densities for U.S.
%	Macroeconomic Variables from a Gaussian Copula Model of Cross-Sectional and Serial Dependence'

%method=2;

switch method
    case 1  % Method 1 using expectations- not accurate.

        n=size(x,1); %input must be an n x 1 vector
        xdash=randsample(x,n);

        E1=mean(abs(x-xobs));
        E2=mean(abs(x-xdash));
        val=E1-E2/2;
        
    case 2  % Method 2 using numerical integration (Appendix B, Smith and Vahey).
         %crpss(t,1)=crps(
         %x=ycdf;
         %xobs=y1;
   
        J=10000;
        step=1/J; a1=step; a2=1-a1;
        a=a1:step:a2;   %   a(1),...,a(J)
     
        x=sort(x);
        xdiff=diff(x);
        eps=mean(xdiff)*0.00001;
        k=find(xdiff<=eps);
        if(~isempty(k)) ;  % jitter observations by eps        
            jitter=normrnd(0,10*eps,size(k));
            x(k)=x(k)+jitter;
        end;
        
        obj=paretotails(x,0,1,'ecdf');  %fit EDF to the data
        q=icdf(obj,a);  %   compute q(1),...,q(J)
        qsval=qs(q,a,xobs); %   compute vector of QS(q(1),xobs),...,QS(q(J),xobs)
        val=mean(qsval);
end

end

