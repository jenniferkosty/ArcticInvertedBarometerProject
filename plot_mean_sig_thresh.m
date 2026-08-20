function plot_mean_sig_thresh(ax, f, dataCell, sigCell, labels, colors, lw, varargin)
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
p.addParameter('Title', ' ')
p.addParameter('FontSize', 12)
p.parse(varargin{:});
S = p.Results;

f = f(:)'; % row
hold(ax,'on');

for i = 1:numel(dataCell)
    x = dataCell{i};
    if isempty(x)
        continue
    end

    %%% Plotting mean
    plot(ax, f, x, 'Color', colors(i,:), 'LineStyle', '-', 'LineWidth', lw, ...
        'DisplayName', string(upper(labels{i})));

    %%% Plotting significance threshold
    if ~isempty(sigCell{i})
        sig = sigCell{i};
        plot(f, sig, 'Color', colors(i,:), 'LineStyle', '--', 'LineWidth', lw, 'HandleVisibility', 'off')
    end
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
fs = S.FontSize;
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
            label_offset_factor = 0.75;   % < 1 shifts label LEFT of the line (in log space), increase toward 1 to move closer, decrease to move further
            text(ax, xpos.*label_offset_factor, yText, txt, 'Rotation', 90, ...
                'HorizontalAlignment','center', 'VerticalAlignment','middle', 'FontSize', fs);
        end
    end
end

if S.ShowLegend
    legend(ax, 'Location', S.LegendLocation);
end

title(S.Title)
box on

set(gca, 'FontSize', fs);

end

