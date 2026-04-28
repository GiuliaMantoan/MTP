function [bQR,bQRbst,bOLS,bOLSbst] = qfe_qr_local_projection_SL(Y_LP,X,quantiles,hz,bst,bstOptions)

% Code to estimate quantile regression
% Can be used to estimate country-specific quantile regression, as well as
% panel quantile regression with quantile-specific country-fixed effects.
% For the latter, need to define 0/1 variables for each country.

% Simon Lloyd and Ed Manuel, March 2022

%% PRELIMINARIES
% =========================================================================
% Get sizes
H  = length(hz);
N  = size(Y_LP,2);          % Number of dependent-variable countries
nQ = length(quantiles);     % Number of quantiles
nX = size(X,2)/N;           % Number of variables

% Bootstrap options
nBoot       = bstOptions.nboot;
blockSize   = bstOptions.blocksize;

% Create empty matrices
bQR     = nan(nX,nQ,H);
bOLS    = nan(nX,H);
bQRbst  = nan(nX,nQ,H,nBoot);
bOLSbst = nan(nX,H,nBoot);

idx = [];
for ii = 1:nX
	idx = [idx ii:nX:N*nX];   
end

%% QR ESTIMATION
% =========================================================================
% Loop over horizons
for h = 1:H
	Ytemp   = rmmissing(Y_LP(:,:,h));       % Removing missing values
    T       = size(Ytemp,1);
    Xestim  = X(1:size(Ytemp,1),:);         % Align dimensions in X
    
    XData = reshape(Xestim(:,idx),T*N,[]);  % Reshape X for estimation
    YData = Ytemp(:);                       % Reshape Y for estimation
    
    % Estimate the QR by quantile
    for qq = 1:nQ
        bQR(:,qq,h) = rq_SL(XData, YData, quantiles(qq));
    % End loop over quantiles
    end
    
    % OLS regression
    bOLS(:,h)       = (XData'*XData)\(XData'*YData);
    
% End loop over horizons    
end

disp('Basic QR Estimation Complete');

%% BOOTSTRAP ESTIMATION
% =========================================================================
if bst == 1
	% Loop over horizons
	for h = 1:H
        % Loop over bootstraps
        for nb = 1:nBoot
            Ytemp   = rmmissing(Y_LP(:,:,h));       % Removing missing values
            T       = size(Ytemp,1);
            Xestim  = X(1:size(Ytemp,1),:);         % Align dimensions in X
            YData = Ytemp(:);                       % Reshape Y for estimation

            draw = overlapping_bst([reshape(YData,T,N) Xestim], blockSize);
            Xdra = draw(:,N+1:end);
            Tbst = size(Xdra,1);
        	Ydra = draw(:,1:N);
            
            Ybst = Ydra(:);
            Xbst = reshape(Xdra(:,idx),Tbst*N,[]);
        
            % Estimate the QR by quantile
            for qq = 1:nQ
                bQRbst(:,qq,h,nb) = rq_SL(Xbst, Ybst, quantiles(qq));
            % End loop over quantiles
            end
            
            % OLS regression
            bOLSbst(:,h,nb)	= (Xbst'*Xbst)\(Xbst'*Ybst);
        % End loop over bootstraps   
        end
    disp(['Bootstraps for Horizon ' num2str(h) ' are Complete']);      
	% End loop over horizons
    end
% End if
end

disp('Boostrap QR Estimation Complete');
% End function
end
