function ax = plot_wavelet(ax, t, x, fs_num, varargin)
% plot_wavelet - Compute and plot normalized wavelet power in a given axes
%
% INPUTS:
%   ax      - axes handle to plot into
%   t       - time vector as datenums
%   x       - time series data vector
%   fs_num  - sampling frequency (cycles/day)
%
% OPTIONAL NAME-VALUE PAIRS:
%   'Title'         - axes title string (default: '')
%   'YLim'          - y-axis limits (default: [min(period_wc) 1050])
%   'YTicks'        - y-axis ticks (default: [0.5 1 2 4 8 16 32 64 128 256 512 1024])
%   'CLim'          - colorbar limits (default: auto)
%   'CbTicks'       - colorbar tick values (default: preset)
%   'CbTickLabels'  - colorbar tick labels (default: preset)
%   'CbLabel'       - colorbar ylabel string (default: 'Normalized Wavelet Power')
%   'FontSize'      - font size (default: 12)
%   'Colormap'      - colormap name or array (default: 'turbo')
%   'ColorScale'    - 'linear' or 'log' (default: 'linear')
%   'ShowColorbar'  - true/false (default: true)

p = inputParser;
p.addParameter('Title',        '',     @(x) ischar(x) || isstring(x));
p.addParameter('YLim',         [],     @isnumeric);  % set after CWT so we know min(period_wc)
p.addParameter('YTicks',       [0.5 1 2 4 8 16 32 64 128 256 512 1024], @isnumeric);
p.addParameter('CLim',         [],     @isnumeric);
p.addParameter('CbTicks',      log10([1/64 1/16 1/4 1 2 4 8 16 64 256 1024]), @isnumeric);
p.addParameter('CbTickLabels', {"1/64","1/16","1/4","1","2","4","8","16","64","256","1024"}, @(x) iscell(x) || isstring(x));
p.addParameter('CbLabel',      {'Normalized', 'Wavelet Power'}, @(x) ischar(x) || isstring(x));
p.addParameter('FontSize',     12,     @isnumeric);
p.addParameter('Colormap',     'turbo');
p.addParameter('ColorScale',   'linear', @(x) ismember(x, {'linear','log'}));
p.addParameter('ShowColorbar', true,   @(x) islogical(x) || isnumeric(x));
p.parse(varargin{:});
S = p.Results;

% Formatting data
x = x(:);
t = t(:);
good   = isfinite(x);
starts = find(diff([0; good]) == 1);
ends   = find(diff([good; 0]) == -1);
ind    = starts(1):ends(end);
x = x(ind);
t = t(ind);
good = isfinite(x);
x = interp1(t(good), x(good), t);

% --- Compute CWT ---
[wt_raw, f_wc, coi] = cwt(x, fs_num);
normPow        = abs(wt_raw).^2 ./ var(x);
period_wc      = 1 ./ f_wc;       % convert frequency to period [days]
coi_period     = 1 ./ coi;        % convert COI frequency to period [days]
T              = yearfrac(datenum(t));      % convert datenums to fractional years

% --- Set YLim default now that we have period_wc ---
if isempty(S.YLim)
    S.YLim = [min(period_wc) 1050];
end

% --- Main pcolor plot ---
hold(ax, 'on');
pcolor(ax, T(:)', period_wc(:), normPow);
shading(ax, 'interp');

% --- Axes formatting ---
set(ax, 'YDir',       'reverse');
set(ax, 'TickDir',    'out');
set(ax, 'YScale',     'log');
set(ax, 'FontSize',   S.FontSize);
set(ax, 'YMinorTick', 'off');
set(ax, 'YTickMode',  'manual');

ylim(ax, S.YLim);
xlim(ax, [T(1) T(end)]);
yticks(ax, S.YTicks);
ylabel(ax, 'Period [days]');

if strlength(S.Title) > 0
    title(ax, S.Title);
end

% --- Colormap ---
colormap(ax, S.Colormap);
ax.ColorScale = S.ColorScale;

if ~isempty(S.CLim)
    clim(ax, S.CLim);
end

% --- Colorbar ---
if S.ShowColorbar
    cb                = colorbar(ax);
    cb.Ticks          = S.CbTicks;
    cb.TickLabels     = S.CbTickLabels;
    cb.TicksMode      = 'manual';
    cb.TickLabelsMode = 'manual';
    ylabel(cb, S.CbLabel, 'Rotation', 270, 'FontSize', S.FontSize);
end

% --- Cone of influence ---
coi_plot = coi_period(:)';
yl       = ylim(ax);
patch(ax, [T(:)', fliplr(T(:)')], ...
    [coi_plot, yl(2) .* ones(1, numel(T))], ...
    'w', 'FaceAlpha', 0.4, 'EdgeColor', 'none');
plot(ax, T(:)', coi_plot, 'w--', 'LineWidth', 1);

end