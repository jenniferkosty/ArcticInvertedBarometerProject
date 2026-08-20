function plot_mean_std_band(ax, f, dataCell, labels, colors, lw, varargin)
% plot_mean_std_band
% dataCell: cell array where dataCell{i} is [nSegments x nFreq]
% labels:   cell array of strings/chars for legend entries

p = inputParser;
p.addParameter('XLog', true, @(x)islogical(x) || isnumeric(x));
p.addParameter('YLog', false, @(x)islogical(x) || isnumeric(x));
p.addParameter('XLim', [min(f) max(f)]);
p.addParameter('YLim', []);
p.addParameter('YTicks', []);
p.addParameter('XLabel', "");
p.addParameter('YLabel', "");
p.addParameter('ShowLegend', false);
p.addParameter('LegendLocation', 'northeast');
p.addParameter('ShowXTicks', true);
p.addParameter('TiePeriods', [], @(x)isnumeric(x));
p.addParameter('TieLabels', {}, @(x)iscell(x) || isstring(x));
p.addParameter('ShowTieText', false);
p.addParameter('isPhase', 0);
p.addParameter('PhaseBandPrct', [16 84], @(x)isnumeric(x) && numel(x)==2); % central band in wrapped deviation

p.parse(varargin{:});
S = p.Results;

f = f(:)'; % row
hold(ax,'on');

for i = 1:numel(dataCell)
    x = dataCell{i};
    if isempty(x)
        continue
    end

    if S.isPhase %%% Taking mean over the unit circle 

        % Circular mean per frequency (radians)
        mu_r = atan2(mean(sin(x),1,'omitnan'), mean(cos(x),1,'omitnan'));  % [-pi, pi]

        % Make mean continuous across frequency to avoid branch-cut artifacts
        mu_r_u = unwrap(mu_r, pi, 2);  % unwrap along frequency

        % Deviations from *wrapped* mean, then wrap deviations to [-pi,pi]
        d = wrapToPi(x - mu_r);        % deviations in [-pi,pi]

        pr = S.PhaseBandPrct;
        dlo = prctile(d, pr(1), 1);    % radians
        dhi = prctile(d, pr(2), 1);    % radians

        % Build band around continuous mean (degrees)
        mu = rad2deg(mu_r_u);
        lo = rad2deg(mu_r_u + dlo);
        hi = rad2deg(mu_r_u + dhi);
        
    else %%% Regular mean and standard deviation
        mu = mean(x, 1, 'omitnan');
        sd = std(x, 0, 1, 'omitnan');
        lo = mu - sd;
        hi = mu + sd;
    end

    %%% Identifying start/stop points
    ok = isfinite(f) & isfinite(lo) & isfinite(hi);
    d = diff([false, ok, false]);
    starts = find(d == 1);
    stops  = find(d == -1) - 1;

    for k = 1:numel(starts)
        idx = starts(k):stops(k);
        shaded_f = [f(idx), fliplr(f(idx))];
        shaded_x = [lo(idx), fliplr(hi(idx))];

        fill(ax, shaded_f, shaded_x, colors(i,:), ...
            'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility','off');
    end

    %%%% Plotting mean
    plot(ax, f, mu, 'Color', colors(i,:), 'LineWidth', lw, ...
        'DisplayName', string(upper(labels{i})));
end

% Axes formatting
xlim(ax, S.XLim);
if ~isempty(S.YLim), ylim(ax, S.YLim); end
if ~isempty(S.YTicks), yticks(ax, S.YTicks); end

if S.XLog, set(ax,'XScale','log'); else, set(ax,'XScale','linear'); end
if S.YLog, set(ax,'YScale','log'); else, set(ax,'YScale','linear'); end

if strlength(S.XLabel) > 0, xlabel(ax, S.XLabel); end
if ~(isempty(S.YLabel) || (isstring(S.YLabel) && strlength(S.YLabel)==0))
    ylabel(ax, S.YLabel);
end

grid(ax,'on');

if ~S.ShowXTicks
    ax.XTickLabel = [];
end

% Tie lines (optionally with text)
if ~isempty(S.TiePeriods)
    xlineY = ylim(ax);
    yText = xlineY(1) + 0.25*(xlineY(2)-xlineY(1));
    for k = 1:numel(S.TiePeriods)
        xpos = 1./S.TiePeriods(k);
        xline(ax, xpos, 'k', 'LineWidth', 1, 'HandleVisibility','off', 'Layer','bottom');
        if S.ShowTieText
            if ~isempty(S.TieLabels)
                txt = S.TieLabels{k};
            else
                txt = sprintf('%g days', S.TiePeriods(k));
            end
            text(ax, xpos/1.1, yText, txt, 'Rotation', 90, ...
                'HorizontalAlignment','center', 'VerticalAlignment','middle');
        end
    end
end

if S.ShowLegend
    legend(ax, 'Location', S.LegendLocation);
end

set(gca, 'FontSize', 12);

end