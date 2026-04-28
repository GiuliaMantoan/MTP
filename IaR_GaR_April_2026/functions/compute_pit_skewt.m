%1. compute_pit_skewt.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [pit, origin_dates, realisation_dates] = compute_pit_skewt( ...
        actual_y, dateNumeric, StartIns, ...
        lc_skt, sc_skt, sh_skt, df_skt, ...
        h, tauGrid)

% COMPUTE_PIT_SKEW_T
%   Probability Integral Transforms (PITs) for a skew-t predictive density.
%
% INPUTS
%   actual_y     : T_full x 1 vector of realised outcomes (e.g. quarterly GDP or inflation)
%   dateNumeric  : T_full x 1 datenums matching actual_y
%   StartIns     : datenum of first in-sample origin used in skew-t arrays
%   lc_skt       : T_ins x H matrix, location parameter for skew-t
%   sc_skt       : T_ins x H matrix, scale
%   sh_skt       : T_ins x H matrix, skewness
%   df_skt       : T_ins x H matrix, degrees of freedom
%   h            : forecast horizon in periods ahead (1..H)
%   tauGrid      : OPTIONAL vector of quantile levels used to approximate
%                  the skew-t CDF (default 0.001:0.001:0.999)
%
% OUTPUTS
%   pit              : N x 1 vector of PITs in [0,1]
%                      N = T_ins - h (last h origins have no realised y_{t+h})
%   origin_dates     : N x 1 datetime of forecast origins (info-set dates)
%   realisation_dates: N x 1 datetime of realisations (dates of y_{t+h})
%
% LOGIC
%   For each origin t in in-sample window:
%       - we have skew-t parameters lc_skt(t,h), sc_skt(t,h), sh_skt(t,h), df_skt(t,h)
%       - we look at the realised y_{t+h}
%       - we compute PIT_t = F_t( y_{t+h} ), approximating F using qskt on a dense tauGrid.

    % ---------- defaults ----------
    if nargin < 10 || isempty(tauGrid)
        tauGrid = (0.001:0.001:0.999)';  % 999 points
    else
        tauGrid = tauGrid(:);
    end

    % ---------- align StartIns with actual_y / dateNumeric ----------
    idx_ins_start = find(dateNumeric >= StartIns, 1, 'first');
    if isempty(idx_ins_start)
        error('StartIns not found or beyond the actual_y sample.');
    end

    T_full = numel(actual_y);
    T_ins  = size(lc_skt,1);       % should match idx_end - idx_ins_start + 1
    H      = size(lc_skt,2);

    if h < 1 || h > H
        error('Requested horizon h=%d is outside [1,%d].', h, H);
    end

    % we assume lc_skt(t,:) corresponds to origin at index:
    %   idx_origin = idx_ins_start + (t-1)  in actual_y / dateNumeric
    idx_end = idx_ins_start + T_ins - 1;
    if idx_end > T_full
        error('Inconsistent sizes: skew-t arrays longer than actual series.');
    end

    % Last origin that has a realised y_{t+h} within the sample:
    max_origin_row = T_ins - h;   % we need origin_row + h <= idx_end
    if max_origin_row < 1
        error('Not enough data to compute PITs for horizon h=%d.', h);
    end

    pit               = NaN(max_origin_row, 1);
    origin_dates      = NaT(max_origin_row, 1);
    realisation_dates = NaT(max_origin_row, 1);

    % ---------- main loop over origins ----------
    for t_row = 1:max_origin_row

        % origin index in actual_y / dateNumeric
        idx_origin = idx_ins_start + (t_row - 1);
        % realisation index: origin + h
        idx_real   = idx_origin + h;

        y_real = actual_y(idx_real);
        if isnan(y_real)
            continue;  % PIT stays NaN
        end

        mu_t  = lc_skt(t_row, h);
        sig_t = sc_skt(t_row, h);
        a_t   = sh_skt(t_row, h);
        nu_t  = df_skt(t_row, h);

        if any(isnan([mu_t,sig_t,a_t,nu_t])) || sig_t <= 0
            continue;
        end

        % ---------- approximate CDF via qskt + interpolation ----------
        q_grid = qskt(tauGrid, mu_t, sig_t, a_t, nu_t);  % same size as tauGrid

        % Ensure q_grid is strictly increasing to avoid issues in interp1
        [q_grid_sorted, sort_idx] = sort(q_grid);
        tau_sorted = tauGrid(sort_idx);

        % approximate PIT as tau such that q(tau) == y_real
        pit_t = interp1(q_grid_sorted, tau_sorted, y_real, 'linear', 'extrap');

        % clip to [0,1] for numerical safety
        pit_t = max(0, min(1, pit_t));

        pit(t_row)               = pit_t;
        origin_dates(t_row)      = datetime(dateNumeric(idx_origin), 'ConvertFrom', 'datenum');
        realisation_dates(t_row) = datetime(dateNumeric(idx_real),   'ConvertFrom', 'datenum');
    end
end


