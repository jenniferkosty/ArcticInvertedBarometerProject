function [ax1, ax2] = plot_wavelet_transfer(ax1, ax2, t, x, y, varargin)
% plot_wavelet_transfer - Compute and plot wavelet transfer function (gain
%                         and phase) between two time series
%
% The wavelet transfer function is defined as:
%   H(s,t) = smooth(Wxy) / smooth(Wxx)
% where Wxy is the cross-wavelet spectrum and Wxx is the power of x.
% Gain = |H|, Phase = angle(H) in degrees.
%
% INPUTS:
%   ax1     - axes handle for gain panel
%   ax2     - axes handle for phase panel
%   t       - time vector as datetimes or datenums
%   x       - input time series (forcing)
%   y       - output time series (response)
%
% OPTIONAL NAME-VALUE PAIRS:
%   'Title'           - title string for gain panel (default: '')
%   'PhaseTitle'      - title string for phase panel (default: 'Transfer Function Phase')
%   'YLim'            - y-axis limits (default: [min(period) 1050])
%   'YTicks'          - y-axis ticks (default: [0.5 1 2 4 8 16 32 64 128 256 512 1024])
%   'CLim'            - colorbar limits for gain panel (default: auto)
%   'FontSize'        - font size (default: 12)
%   'Colormap'        - colormap for gain panel (default: 'turbo')
%   'SmoothWidth'     - Gaussian smoothing width in time steps (default: 10)
%   'RunMonteCarlo'   - true = run MC significance testing (default: false)
%   'NSurrogates'     - number of surrogates (default: 200)
%   'SigLevel'        - significance level as percentile (default: 95)
%   'SigMask'         - precomputed logical [nScales x nTime] mask
%   'ShowSigContour'  - overlay significance contour (default: true)
%   'LogGain'         - plot gain on log10 scale (default: false)

p = inputParser;
p.addParameter('Title',           '',       @(x) ischar(x) || isstring(x));
p.addParameter('PhaseTitle',      'Transfer Function Phase', @(x) ischar(x) || isstring(x));
p.addParameter('PlotPhase',       true,     @(x) islogical(x) || isnumeric(x));
p.addParameter('YLim',            [],       @isnumeric);
p.addParameter('YTicks',          [0.5 1 2 4 8 16 32 64 128 256 512 1024], @isnumeric);
p.addParameter('CLim',            [],       @isnumeric);
p.addParameter('FontSize',        12,       @isnumeric);
p.addParameter('Colormap',        'turbo',  @(x) ischar(x) || isstring(x) || isnumeric(x));
p.addParameter('SmoothWidth',     10,       @isnumeric);
p.addParameter('RunMonteCarlo',   false,    @(x) islogical(x) || isnumeric(x));
p.addParameter('NSurrogates',     200,      @isnumeric);
p.addParameter('SigLevel',        95,       @isnumeric);
p.addParameter('SigMask',         [],       @(x) islogical(x) || isnumeric(x));
p.addParameter('ShowSigContour',  true,     @(x) islogical(x) || isnumeric(x));
p.addParameter('LogGain',         false,    @(x) islogical(x) || isnumeric(x));
p.addParameter('MaxPeriod',       500,     @isnumeric);  
p.parse(varargin{:});
S = p.Results;

% =========================================================================
% --- Prepare data ---
% =========================================================================
t = t(:);
x = x(:);
y = y(:);

%%% Remove NaNs from x
good   = isfinite(x) & isfinite(y);
starts = find(diff([0; good]) == 1);
ends   = find(diff([good; 0]) == -1);
ind    = starts(1):ends(end);
x      = x(ind);
t      = t(ind);
good   = isfinite(x);
x      = interp1(t(good), x(good), t);

%%% Remove NaNs from y
y      = y(ind);
good = isfinite(y);
y      = interp1(t(good), y(good), t);

%%% Sampling frequency
t_datenum = datenum(t);
dt        = t_datenum(2) - t_datenum(1);
fs_num    = 1 / dt;

%%% Time axis in fractional years
T = yearfrac(t_datenum);

% =========================================================================
% --- Compute auto and cross-spectra
% =========================================================================
[~, wcs, f_wc, coi] = wcoherence(x, y, fs_num);
[~, Wxx, ~, ~] = wcoherence(x, x, fs_num);

% =========================================================================
% --- Compute transfer function ---
% =========================================================================

H_gain      = abs(wcs) ./ Wxx;        % Pxy / Pxx, same normalisation throughout
H_phase = rad2deg(angle(wcs));

period = 1 ./ f_wc;   % [nScales x 1]

coi_period = 1 ./ coi;

%%% Set YLim
if isempty(S.YLim)
    S.YLim = [min(period) S.MaxPeriod];
end

if S.LogGain
    H_gain = log10(H_gain);
end

% =========================================================================
% --- Monte Carlo significance testing ---
% =========================================================================
sig_mask = S.SigMask;

if S.RunMonteCarlo && isempty(sig_mask)

    N       = numel(x);
    nScales = numel(fx);
    nSur    = S.NSurrogates;

    fprintf('Running Monte Carlo with %d surrogates...\n', nSur);

    x0 = x - mean(x, 'omitnan');
    y0 = y - mean(y, 'omitnan');

    gain_sur = zeros(nScales, nSur);

    for i = 1:nSur
        if mod(i, 50) == 0
            fprintf('  Surrogate %d / %d\n', i, nSur);
        end

        %%% Phase-randomize y only
        Y           = fft(y0);
        Nhalf       = floor(N/2);
        rand_phases = exp(1i * 2*pi * rand(N, 1));
        rand_phases(1) = 1;

        if mod(N, 2) == 0
            rand_phases(Nhalf+1) = real(rand_phases(Nhalf+1));
        end

        mirror_idx = Nhalf+2 : N;
        if mod(N, 2) == 0
            source_idx = Nhalf : -1 : 2;
        else
            source_idx = Nhalf+1 : -1 : 2;
        end
        rand_phases(mirror_idx) = conj(rand_phases(source_idx));
        yS = real(ifft(Y .* rand_phases));

        %%% CWT of surrogate
        [Wx_s, ~] = cwt(x0, fs_num);
        [Wy_s, ~] = cwt(yS, fs_num);

        if size(Wx_s, 1) == numel(T)
            Wx_s = Wx_s.';
            Wy_s = Wy_s.';
        end

        %%% Surrogate transfer function gain
        Wxy_s    = Wx .* conj(Wy_s);
        Wxx_s    = abs(Wx).^2;

        H_s      = Wxy_s ./ Wxx_s;
        gain_s   = abs(H_s);

        %%% Masking same threshold
        log_Wxx_s = log10(Wxx_s);
        low_power_mask = log_Wxx_s < thresh;

        gain_s(low_power_mask)  = NaN;

        if S.LogGain
            gain_s = log10(gain_s);
        end

        gain_sur(:, i) = prctile(gain_s, S.SigLevel, 2);
    end

    sig_thresh_mc = prctile(gain_sur, S.SigLevel, 2);   % [nScales x 1]
    sig_mask      = H_gain >= sig_thresh_mc;
    fprintf('Monte Carlo complete.\n');

end

% =========================================================================
% --- Panel 1: Gain ---
% =========================================================================
hold(ax1, 'on');

pcolor(ax1, T(:)', period(:), H_gain);
shading(ax1, 'interp');

set(ax1, 'YDir',       'reverse');
set(ax1, 'TickDir',    'out');
set(ax1, 'YScale',     'log');
set(ax1, 'FontSize',   S.FontSize);
set(ax1, 'YMinorTick', 'off');
set(ax1, 'YTickMode',  'manual');

ylim(ax1, S.YLim);
xlim(ax1, [T(1) T(end)]);
yticks(ax1, S.YTicks);
ylabel(ax1, 'Period [days]');

if strlength(string(S.Title)) > 0
    title(ax1, S.Title);
end

colormap(ax1, S.Colormap);
cb1 = colorbar(ax1);
if ~isempty(S.CLim)
    clim(ax1, S.CLim);
end
if S.LogGain
    cb1.Label.String = 'Gain [log_{10}]';
else
    cb1.Label.String = 'Gain |H(s,t)|';
end
cb1.Label.Rotation = 270;
cb1.Label.FontSize = S.FontSize;

%%% COI
coi_plot = coi_period(:)';
yl       = ylim(ax1);
patch(ax1, [T(:)', fliplr(T(:)')], ...
    [coi_plot, yl(2) .* ones(1, numel(T))], ...
    'w', 'FaceAlpha', 0.4, 'EdgeColor', 'none');
plot(ax1, T(:)', coi_plot, 'w--', 'LineWidth', 1);

%%% Significance contour
if S.ShowSigContour && ~isempty(sig_mask)
    contour(ax1, T(:)', period(:), double(sig_mask), [0.5 0.5], 'k', 'LineWidth', 1.5);
end

% =========================================================================
% --- Panel 2: Phase ---
% =========================================================================
if S.PlotPhase

H_phase_masked = H_phase;
if ~isempty(sig_mask)
    H_phase_masked(~sig_mask) = NaN;
end

hold(ax2, 'on');
pcolor(ax2, T(:)', period(:), H_phase_masked);
shading(ax2, 'interp');

set(ax2, 'YDir',       'reverse');
set(ax2, 'TickDir',    'out');
set(ax2, 'YScale',     'log');
set(ax2, 'FontSize',   S.FontSize);
set(ax2, 'YMinorTick', 'off');
set(ax2, 'YTickMode',  'manual');

ylim(ax2, S.YLim);
xlim(ax2, [T(1) T(end)]);
yticks(ax2, S.YTicks);
ylabel(ax2, 'Period [days]');
xlabel(ax2, 'Year');
title(ax2, S.PhaseTitle);

colormap(ax2, hsv);
cb2               = colorbar(ax2);
clim(ax2,         [-180 180]);
cb2.Ticks         = -180:60:180;
cb2.Label.String  = 'Phase [deg]';
cb2.Label.Rotation = 270;
cb2.Label.FontSize = S.FontSize;

%%% COI
coi_plot = coi_period(:)';
yl       = ylim(ax2);
patch(ax2, [T(:)', fliplr(T(:)')], ...
    [coi_plot, yl(2) .* ones(1, numel(T))], ...
    'w', 'FaceAlpha', 0.4, 'EdgeColor', 'none');
plot(ax2, T(:)', coi_plot, 'w--', 'LineWidth', 1);

%%% Significance contour
if S.ShowSigContour && ~isempty(sig_mask)
    contour(ax2, T(:)', period(:), double(sig_mask), [0.5 0.5], 'k', 'LineWidth', 0.5);
end

end

end