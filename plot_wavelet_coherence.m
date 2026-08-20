function [ax1, ax2] = plot_wavelet_coherence(ax1, ax2, t, x, y, varargin)
% plot_wavelet_coherence - Compute and plot wavelet coherence and relative
%                          phase in two given axes
%
% INPUTS:
%   ax1     - axes handle for coherence/shared power panel
%   ax2     - axes handle for relative phase panel (can be [] if PlotPhase = false)
%   t       - time vector as datenums
%   x       - first time series
%   y       - second time series
%
% OPTIONAL NAME-VALUE PAIRS:
%   'fs'              - sampling frequency in cycles/day (default: 1/diff(t(1:2)))
%   'PlotCoherence'   - true = plot coherence, false = plot shared power (default: false)
%   'PlotPhase'       - true = plot relative phase in ax2 (default: true)
%   'Title'           - title string for top panel (default: '')
%   'PhaseTitle'      - title string for phase panel (default: 'Relative Phase (+: x leads y)')
%   'VarName'         - variable name string for titles (default: 'x')
%   'YLim'            - y-axis limits (default: [min(period_wc) 1050])
%   'YTicks'          - y-axis ticks (default: [0.5 1 2 4 8 16 32 64 128 256 512 1024])
%   'CLim'            - colorbar limits for top panel (default: auto)
%   'FontSize'        - font size (default: 12)
%   'Colormap'        - colormap for top panel (default: 'turbo')
%   'RunMonteCarlo'   - true = run MC significance testing (default: false)
%   'NSurrogates'     - number of surrogates for MC test (default: 200)
%   'SigLevel'        - significance level as percentile (default: 95)
%   'SigMask'         - precomputed logical [nScales x nTime] significance mask
%                       (overrides RunMonteCarlo if provided)
%   'ShowSigContour'  - overlay significance contour if SigMask provided (default: true)

p = inputParser;
p.addParameter('fs',             [],       @isnumeric);
p.addParameter('PlotCoherence',  false,    @(x) islogical(x) || isnumeric(x));
p.addParameter('PlotPhase',      true,     @(x) islogical(x) || isnumeric(x));
p.addParameter('Title',          '',       @(x) ischar(x) || isstring(x));
p.addParameter('PhaseTitle',     'Relative Phase (+: x leads y)', @(x) ischar(x) || isstring(x));
p.addParameter('VarName',        'x',      @(x) ischar(x) || isstring(x));
p.addParameter('YLim',           [],       @isnumeric);
p.addParameter('YTicks',         [0.5 1 2 4 8 16 32 64 128 256 512 1024], @isnumeric);
p.addParameter('CLim',           [],       @isnumeric);
p.addParameter('FontSize',       12,       @isnumeric);
p.addParameter('Colormap',       'turbo',  @(x) ischar(x) || isstring(x) || isnumeric(x));
p.addParameter('RunMonteCarlo',  false,    @(x) islogical(x) || isnumeric(x));
p.addParameter('NSurrogates',    200,      @isnumeric);
p.addParameter('SigLevel',       95,       @isnumeric);
p.addParameter('SigMask',        [],       @(x) islogical(x) || isnumeric(x));
p.addParameter('ShowSigContour', true,     @(x) islogical(x) || isnumeric(x));
p.parse(varargin{:});
S = p.Results;

% Formatting data
t = t(:);
x = x(:);
y = y(:);

%%% Remove NaNs from x
good   = isfinite(x) & isfinite(y);
starts = find(diff([0; good]) == 1);
ends   = find(diff([good; 0]) == -1);
ind    = starts(1):ends(end);
x = x(ind);
t = t(ind);
good = isfinite(x);
x = interp1(t(good), x(good), t);

%%% Remove NaNs from y
y = y(ind);
good = isfinite(y);
y = interp1(t(good), y(good), t);

%%% Removing mean
x = x - mean(x, 'omitnan');
y = y - mean(y, 'omitnan');

%%% Compute sampling frequency
t_datenum = datenum(t);
dt        = t_datenum(2) - t_datenum(1);
fs_num    = 1 / dt;

%%% Compute coherence
[wcoh, wcs, f_wc, coi] = wcoherence(x, y, fs_num);
period_wc  = 1 ./ f_wc;
coi_period = 1 ./ coi;
T          = yearfrac(t_datenum);
Phase_rel  = rad2deg(angle(wcs));

%%% Ensure correct orientation
if size(wcoh, 1) == numel(T)
    wcoh      = wcoh.';
    wcs       = wcs.';
    Phase_rel = Phase_rel.';
end

%%% Set ylim
if isempty(S.YLim)
    S.YLim = [min(period_wc) 1050];
end

% =========================================================================
% --- Monte Carlo significance testing ---
% =========================================================================
sig_mask = S.SigMask;   % use precomputed mask if provided

if S.RunMonteCarlo && isempty(sig_mask)
    nSur    = S.NSurrogates;
    N       = numel(x);
    nScales = numel(f_wc);

    fprintf('Running Monte Carlo with %d surrogates...\n', nSur);

    sp_real = abs(wcs) ./ (std(x) .* std(y));
    sp_sur  = zeros(nScales, nSur);

    for i = 1:nSur
        if mod(i, 50) == 0
            fprintf('  Surrogate %d / %d\n', i, nSur);
        end

        Y           = fft(y);
        Nhalf       = floor(N/2);
        rand_phases = exp(sqrt(-1) * 2 * pi * rand(N, 1));

        rand_phases(1) = 1;

        if mod(N, 2) == 0
            rand_phases(Nhalf+1) = real(rand_phases(Nhalf+1));
        end
        mirror_idx = Nhalf+2 : N;

        if mod(N,2) == 0
            source_idx = Nhalf : -1 : 2;
        else
            source_idx = Nhalf+1 : -1 : 2;
        end
        rand_phases(mirror_idx) = conj(rand_phases(source_idx));
        yS = real(ifft(Y .* rand_phases));

        [~, wcs_s, ~, ~] = wcoherence(x, yS, fs_num);

        if size(wcs_s, 1) == numel(T)
            wcs_s = wcs_s.';
        end

        sp_s = abs(wcs_s) ./ (std(x) .* std(y));
        sp_sur(:, i) = prctile(sp_s, S.SigLevel, 2);
    end

    sig_thresh_mc = prctile(sp_sur, S.SigLevel, 2);
    sig_mask      = sp_real >= sig_thresh_mc;
    fprintf('Monte Carlo complete.\n');
end

% =========================================================================
% --- Panel 1: Coherence or Shared Power ---
% =========================================================================
hold(ax1, 'on');

if S.PlotCoherence
    plotData = wcoh;
else
    plotData = abs(wcs) ./ (std(x) .* std(y));
end

pcolor(ax1, T(:)', period_wc(:), plotData);
shading(ax1, 'interp');

set(ax1, 'YDir',       'reverse');
set(ax1, 'TickDir',    'out');
set(ax1, 'YScale',     'log');
set(ax1, 'FontSize',   S.FontSize);
set(ax1, 'YMinorTick', 'off');
set(ax1, 'YTickMode',  'manual');
xlabel(ax1, 'Year');

ylim(ax1, S.YLim);
xlim(ax1, [T(1) T(end)]);
yticks(ax1, S.YTicks);
ylabel(ax1, 'Period [days]');

if strlength(string(S.Title)) > 0
    title(ax1, S.Title);
end

colormap(ax1, S.Colormap);

if S.PlotCoherence
    cb1                = colorbar(ax1);
    clim(ax1,          [0 1]);
    cb1.Ticks          = 0:0.2:1;
    cb1.Label.String   = 'Magnitude Squared Coherence';
    cb1.Label.Rotation = 270;
    cb1.Label.FontSize = S.FontSize;
else
    cb1 = colorbar(ax1);
    if ~isempty(S.CLim)
        clim(ax1, S.CLim);
    end
    cb1.Label.String   = {'Normalized', 'Shared Power'};
    cb1.Label.Rotation = 270;
    cb1.Label.FontSize = S.FontSize;
end

% --- COI panel 1 ---
coi_plot = coi_period(:)';
yl       = ylim(ax1);
patch(ax1, [T(:)', fliplr(T(:)')], ...
    [coi_plot, yl(2) .* ones(1, numel(T))], ...
    'w', 'FaceAlpha', 0.4, 'EdgeColor', 'none');
plot(ax1, T(:)', coi_plot, 'w--', 'LineWidth', 1);

% --- Significance contour panel 1 ---
if S.ShowSigContour && ~isempty(sig_mask)
    contour(ax1, T(:)', period_wc(:), double(sig_mask), [0.5 0.5], 'k', 'LineWidth', 1.5);
end

% =========================================================================
% --- Panel 2: Relative Phase ---
% =========================================================================
if S.PlotPhase

    Phase_masked = Phase_rel;
    if ~isempty(sig_mask)
        Phase_masked(~sig_mask) = NaN;
    end

    hold(ax2, 'on');
    pcolor(ax2, T(:)', period_wc(:), Phase_masked);
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
    cb2                = colorbar(ax2);
    clim(ax2,          [-180 180]);
    cb2.Ticks          = -180:60:180;
    cb2.Label.String   = 'Phase [deg]';
    cb2.Label.Rotation = 270;
    cb2.Label.FontSize = S.FontSize;

    % --- COI panel 2 ---
    coi_plot = coi_period(:)';
    yl       = ylim(ax2);
    patch(ax2, [T(:)', fliplr(T(:)')], ...
        [coi_plot, yl(2) .* ones(1, numel(T))], ...
        'w', 'FaceAlpha', 0.4, 'EdgeColor', 'none');
    plot(ax2, T(:)', coi_plot, 'w--', 'LineWidth', 1);

    % --- Significance contour panel 2 ---
    if S.ShowSigContour && ~isempty(sig_mask)
        contour(ax2, T(:)', period_wc(:), double(sig_mask), [0.5 0.5], 'k', 'LineWidth', 0.5);
    end

end

end