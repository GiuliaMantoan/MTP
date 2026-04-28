function [fittedParams, fittedQuantiles, PDFx, fittedPDFy, fittedMoments,fittedMomentsPM] = fit_skewt_to_quantiles_all(AllQuantiles,quantiles,quantilesToFit,countriesToFit,horsToFit,maxdf)

% This version takes as inputs quantiles estimated across time, countries
% and horizons

% fits skew-t distribution to quantiles (see
% https://www.newyorkfed.org/medialibrary/media/research/staff_reports/sr794.pdf?la=en,
% p.9 and http://azzalini.stat.unipd.it/SN/ for the matlab code and other
% resources. This function allows both exactly-identified and
% over-identified fitting.

if length(quantilesToFit)<4
    
    error('Not enough quantiles to fit parameters')
    
end

% Preliminaries
periods     = size(AllQuantiles,1);
nQ          = length(quantilesToFit);
npoints     = 400; % how many points to fit PDF to
nctry       = length(countriesToFit);
nhor        = length(horsToFit); 
maxdf       = 30;%maxdf; % what is max for degrees of freedom paramteter 

quantileIdxToFit = arrayfun(@(q) find(abs(quantiles - q) < 1e-8, 1),quantilesToFit);
fittedParams        = zeros(periods,4,nctry,nhor);
fittedQuantiles     = zeros(periods,nQ,nctry,nhor);
fittedPDFy          = zeros(periods,npoints,nctry,nhor);
fittedMoments       = zeros(periods,4,nctry,nhor);

% Plagborg-Moller Bounds
lb = [     -20,     0,   -30];
ub = [      20,    50,   30];

%     Adrian (2019) bounds
% lb = [     -20,     0,   -30];
% ub = [      20,    50,    30];

% Loop over countries
for ctry=1:length(countriesToFit)
    % Loop of horizons
    for hor=1:length(horsToFit)
        % Loop over periods
        for t = 1:periods
            % Fitted Quantiles at the period in question
            lpQuantiles=squeeze(AllQuantiles(t,countriesToFit(ctry),:,horsToFit(hor)));
    
            % Set initial conditions to align with PlagMol
            iqn = norminv(0.75) - norminv(0.25);    
            med = find( abs(quantiles - 0.5) <  1e-8 );
            lc0 = AllQuantiles(t,countriesToFit(ctry),med,horsToFit(hor)); %median         
            uq  = find(quantiles==0.75);
            lq  = find(quantiles==0.25);
            sc0 = (lpQuantiles(uq) - lpQuantiles(lq)) / iqn;  % ratio of IQR vs normal IQR     
            sh0 = 0;
            X0 = [lc0, sc0, sh0];
            
            %% Solve for the parameters: lsqnonlin
            par = NaN(maxdf,3);
            ssq = NaN(maxdf,1);
            for df = 5:maxdf % need df>4, hence start from 4; trial and error proves don't need to go above 30 (also true in Adrian and PM)
                [par(df,:), ssq(df)] = lsqnonlin( @(x) ((lpQuantiles(quantileIdxToFit)' - qskt(quantiles(quantileIdxToFit),x(1),x(2),x(3),df))),...
                    X0(1:3),lb,ub,  optimoptions('lsqnonlin','Display','off')); % optimoptions added by SM aug 25
                %[par1(df,:), ssq1(df)] = fmincon( @(x) sum((lpQuantiles(quantileIdxToFit)' - qskt(quantiles(quantileIdxToFit),x(1),x(2),x(3),df)).^2),...
                %    X0(1:3),[],[],[],[],lb,ub);
            end
            X = NaN(1,4);
            [~,X(4)] = min(ssq);
            X(1:3) = par(X(4),:);
            fittedParams(t,:,ctry,hor) = X;
            
            %% Remaining Quantites
            fittedQuantiles(t,:,ctry,hor) = qskt(quantiles(quantileIdxToFit),fittedParams(t,1,ctry,hor),fittedParams(t,2,ctry,hor),fittedParams(t,3,ctry,hor),fittedParams(t,4,ctry,hor));
    
            fittedCumulants(t,:,ctry,hor) = skt_cumulants(fittedParams(t,1,ctry,hor),fittedParams(t,2,ctry,hor), fittedParams(t,3,ctry,hor), fittedParams(t,4,ctry,hor));

            rangeToFit = linspace(min(lpQuantiles(:))-6,max(lpQuantiles(:))+6,npoints);
            deltaYY    = (rangeToFit(end)-rangeToFit(1))/npoints;

            fittedPDFy(t,:,ctry,hor) = dskt(rangeToFit,fittedParams(t,1,ctry,hor),fittedParams(t,2,ctry,hor),fittedParams(t,3,ctry,hor),fittedParams(t,4,ctry,hor));
            PDFx(t,:,ctry,hor)= rangeToFit;
            
            % Previous approach to calculating moments using the cumulants

            fittedMoments(t,1,ctry,hor)=fittedCumulants(t,1,ctry,hor);
            fittedMoments(t,2,ctry,hor)=fittedCumulants(t,2,ctry,hor);
            fittedMoments(t,3,ctry,hor)=fittedCumulants(t,3,ctry,hor)/(fittedCumulants(t,2,ctry,hor).^(3/2));
            fittedMoments(t,4,ctry,hor)=fittedCumulants(t,4,ctry,hor)./(fittedCumulants(t,2,ctry,hor).^2)+3;
            
            % Alternate approach to fitting moments (Plagborg-Moller,
            % 2021):
            
            fittedMomentsPM (t,1,ctry,hor) = sum(PDFx(t,:,ctry,hor).*fittedPDFy(t,:,ctry,hor)*deltaYY);
            fittedMomentsPM (t,2,ctry,hor) = sum(((PDFx(t,:,ctry,hor)-  fittedMomentsPM (t,1,ctry,hor)).^2).*fittedPDFy(t,:,ctry,hor)*deltaYY);
            fittedMomentsPM (t,3,ctry,hor) = sum(((PDFx(t,:,ctry,hor)-  fittedMomentsPM (t,1,ctry,hor)).^3).*fittedPDFy(t,:,ctry,hor)*deltaYY)/(fittedMomentsPM (t,2,ctry,hor)^(3/2));
            fittedMomentsPM (t,4,ctry,hor) = sum(((PDFx(t,:,ctry,hor)-  fittedMomentsPM (t,1,ctry,hor)).^4).*fittedPDFy(t,:,ctry,hor)*deltaYY)/(fittedMomentsPM (t,2,ctry,hor)^(2));


        end
    end
end

end