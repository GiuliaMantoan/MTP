function [bQR,bQRbst,bOLS,bOLSbst] = qfe_qr_local_projection_SL_final(Y_LP,X,quantiles,hz,bst,bstOptions, covid, exo)

% Code to estimate quantile regression
% Can be used to estimate country-specific quantile regression, as well as
% panel quantile regression with quantile-specific country-fixed effects.
% For the latter, need to define 0/1 variables for each country.

% Simon Lloyd and Ed Manuel, March 2022
% modified by Simone Maso to account for COVID dummy

%% PRELIMINARIES
% =========================================================================

% add covid dummy to the vector of explanatory variables
if covid == 1
    X = [X, exo];
end

% Get sizes
H  = length(hz);
N  = size(Y_LP,2);          % Number of dependent-variable countries
nQ = length(quantiles);     % Number of quantiles
nX = size(X,2)/N;           % Number of variables

% Bootstrap options
nBoot       = bstOptions.nboot;
blockSize   = bstOptions.blocksize;

% Create empty matrices
bQR     = zeros(nX,nQ,H);   
bOLS    = zeros(nX,H);     
bQRbst  = zeros(nX,nQ,H,nBoot);    
bOLSbst = zeros(nX,H,nBoot);    

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

    %% Modify by SM June 2025
    % we need to account for the estimation time close to the covid period 
    % for which - at certain horizons, the vector of exogenous variables
    % will be collinear with the constant
    
    final_idx = size(idx,2); % size of final index
    flag = 0; % 0 normal size | 1 does not account fro exo var

    if all( Xestim(:,end) == 0 ) 
        Xestim = Xestim(:,1:end-1);
        flag = 1;
        final_idx = size(idx,2) -1;
    end 
    
    XData = reshape(Xestim(:,idx(1:final_idx)),T*N,[]);  % Reshape X for estimation
    YData = Ytemp(:);                       % Reshape Y for estimation
    
    % Estimate the QR by quantile
    for qq = 1:nQ

    %% Modify by SM June 2025

        if flag == 1
            bQR(1:end-1,qq,h) = rq_SL(XData, YData, quantiles(qq));
        else
            bQR(:,qq,h) = rq_SL(XData, YData, quantiles(qq));
        end

    % End loop over quantiles
    end
    
    % OLS regression
    if flag == 1
        bOLS(1:end-1,h)  = (XData'*XData)\(XData'*YData);
    else
        bOLS(:,h)       = (XData'*XData)\(XData'*YData);
    end

% End loop over horizons    
end

% disp('Basic QR Estimation Complete');

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

            final_idx = size(idx,2); % size of final index
            flag = 0; % 0 normal size | 1 does not account fro exo var
            if all( Xestim(:,end) == 0 )
                Xestim = Xestim(:,1:end-1);
                flag = 1;
                final_idx = size(idx,2) -1;
            end

            draw = overlapping_bst_SL_exo_var([reshape(YData,T,N) Xestim], blockSize);
            Xdra = draw(:,N+1:end);
            Tbst = size(Xdra,1);
        	Ydra = draw(:,1:N);
            
            Ybst = Ydra(:);
            Xbst = reshape(Xdra(:,idx(1:final_idx)),Tbst*N,[]);
        
            % Estimate the QR by quantile
            for qq = 1:nQ
                if flag == 1
                    bQRbst(1:end-1,qq,h,nb) = rq_SL(Xbst, Ybst, quantiles(qq));
                else
                    bQRbst(:,qq,h,nb) = rq_SL(Xbst, Ybst, quantiles(qq));
                end
            % End loop over quantiles
            end
            
            % OLS regression
            if flag == 1
                bOLSbst(1:end-1,h,nb)	= (Xbst'*Xbst)\(Xbst'*Ybst);
            else 
                bOLSbst(:,h,nb)	=         (Xbst'*Xbst)\(Xbst'*Ybst);
            end
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
