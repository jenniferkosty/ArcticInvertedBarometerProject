%% Script to visualize spectral estimates for OBP, SLP

clc; clear;

%% Loading mooring data

moorings_all = {'a', 'b', 'd'};
dt_string = '30min';
for i = 1:length(moorings_all)

    %%% Selecting mooring
    mooring_no = moorings_all{i};

    %%% Loading BPR data
    data = load('C:\Users\jak279\OneDrive - Yale University\Research\Data\OBP_interp\BPR_mooring_' + string(mooring_no) + '_' + string(dt_string) + '.mat'); % BPR
    bpr.datetime = data.bpr.datetime;
    bpr_tmp = data.bpr.detided;
    %%% Converting to units of m
    if ~strcmp(mooring_no, 'np') %%% NP data is already in m
        bpr_tmp = -gsw_z_from_p(bpr_tmp, data.bpr.lat);
    end
    bpr.(string(mooring_no) + '_detided') = bpr_tmp;

    bpr_tmp = data.bpr.pressure_anomaly;
    %%% Converting to units of m
    if ~strcmp(mooring_no, 'np') %%% NP data is already in m
        bpr_tmp = -gsw_z_from_p(bpr_tmp, data.bpr.lat);
    end
    bpr.(string(mooring_no) + '_pressure_anomaly') = bpr_tmp;
    bpr.(string(mooring_no) + '_start') = data.bpr.start;
    bpr.(string(mooring_no) + '_end') = data.bpr.end;
    if strcmp(mooring_no, 'np') == 1
        bpr.(string(mooring_no) + '_duration') = data.bpr.duration;
    end

    %%% Loading SLP data
    data = load('C:\Users\jak279\OneDrive - Yale University\Research\Data\OBP_interp\SLP_mooring_' + string(mooring_no) + '_' + string('6hr') + '.mat'); % SLP
    slp.datetime = data.slp.datetime;
    slp.(string(mooring_no)) = NaN(size(data.slp.pressure));
    slp_tmp = data.slp.pressure;
    slp_tmp = -gsw_z_from_p(slp_tmp, data.slp.actual_lat);
    % slp_tmp = slp_tmp * 1e4; % Converting to Pa
    % slp_tmp = slp_tmp ./ (rho0 * g); % Converting to m
    for u = 1:length(bpr.(string(mooring_no) + '_start'))
        t1 = bpr.(string(mooring_no) + '_start')(u);
        t2 = bpr.(string(mooring_no) + '_end')(u);
        ind = data.slp.datetime >= t1 & data.slp.datetime <= t2;
        slp.(string(mooring_no))(ind) = slp_tmp(ind) - mean(slp_tmp(ind), 'omitnan');
    end

end

%% Computing average spectral estimates

fs = 12;
fig = figure('Position', [10 10 1000 500]);
tiledlayout(1,2, 'TileSpacing', 'tight', 'Padding', 'tight')
colors = orderedcolors('gem');
lw = 2;

%%% Level of smoothing
P = 4;
numSegs = 8;
use_pwelch = 1;

%%%%%%%%%%%%%%%%%%%%
%%% Raw OBP data %%%
%%%%%%%%%%%%%%%%%%%%

nexttile()

%%% Computing spectral estimate for each deployment
clear sxx_all
for i = 1:length(moorings_all)
    mooring_no = moorings_all{i};
    sxx_all{i} = NaN(length(bpr.(string(mooring_no) + '_start')), length(f));

    for u = 1:length(bpr.(string(mooring_no) + '_start'))

        %%% Isolating annual segments
        t1 = bpr.(string(mooring_no) + '_start')(u);
        t2 = bpr.(string(mooring_no) + '_end')(u);
        idx = bpr.datetime >= t1 & bpr.datetime <= t2;
        t = bpr.datetime(idx);
        x = bpr.(string(mooring_no) + '_pressure_anomaly')(idx);

        %%% Skipping if segment is empty
        if isempty(x)
            sxx_all{i}(u,:) = NaN(size(f));
            continue
        end
        
        %%% Computing spectra
        if use_pwelch == 1
            [f_tmp, sxx_tmp] = spectral_estimate_pwelch(t, x, 'NumSegments', numSegs, 'MakeFigure', 0);
        else
            [f_tmp, sxx_tmp] = spectral_estimate_jlab(t, x, 'P', P, 'MakeFigure', 0);
        end

        %%% Skipping if segment wasn't long enough
        if isnan(f_tmp)
            sxx_all{i}(u,:) = NaN(size(f));
            continue
        end

        %%% Interpolating to common frequency grid for averaging
        idx = ~isnan(sxx_tmp);
        sxx_all{i}(u,:) = interp1(f_tmp(idx), sxx_tmp(idx), f, 'linear', NaN);

    end
    
end

%%% Plotting spectral estimate
hold on
f = f(:)';

for i = 1:length(moorings_all)
    mooring_no = moorings_all{i};
    x = sxx_all{i};

    mu = mean(log10(x), 1, 'omitnan');
    sd = std(log10(x), 0, 1, 'omitnan');

    lo = 10.^(mu - sd);
    hi = 10.^(mu + sd);

    % valid where both bounds are positive and finite (good for ylog)
    ok = isfinite(f) & isfinite(lo) & isfinite(hi) & (lo > 0) & (hi > 0);

    % Find contiguous valid segments
    d = diff([false, ok, false]);
    starts = find(d == 1);
    stops  = find(d == -1) - 1;

    % Draw one fill per segment (so it stops/restarts)
    for k = 1:numel(starts)
        idx = starts(k):stops(k);
        shaded_f = [f(idx), fliplr(f(idx))];
        shaded_x = [lo(idx), fliplr(hi(idx))];

        fill(shaded_f, shaded_x, colors(i,:), ...
            'FaceAlpha', 0.3, ...
            'EdgeColor', 'none', 'HandleVisibility', 'off');
    end

    plot(f, 10.^mu, 'LineStyle', '--', 'Color', colors(i,:), 'LineWidth', lw, ...
        'HandleVisibility', 'off');
end

xlog; ylog;
xlim([1/180 max(f)])
xlabel('Frequency [cycles/day]')
ylabel('Power Spectral Density')
grid on

%%%%%%%%%%%%%%%%%%%%%%%%
%%% Detided OBP data %%%
%%%%%%%%%%%%%%%%%%%%%%%%

%%% Computing spectral estimate for each deployment
clear sxx_all
for i = 1:length(moorings_all)
    mooring_no = moorings_all{i};
    sxx_all{i} = NaN(length(bpr.(string(mooring_no) + '_start')), length(f));

    for u = 1:length(bpr.(string(mooring_no) + '_start'))

        %%% Isolating annual segments
        t1 = bpr.(string(mooring_no) + '_start')(u);
        t2 = bpr.(string(mooring_no) + '_end')(u);
        idx = bpr.datetime >= t1 & bpr.datetime <= t2;
        t = bpr.datetime(idx);
        x = bpr.(string(mooring_no) + '_detided')(idx);

        %%% Skipping if segment is empty
        if isempty(x)
            sxx_all{i}(u,:) = NaN(size(f));
            continue
        end
        
        %%% Computing spectra
        if use_pwelch == 1
            [f_tmp, sxx_tmp] = spectral_estimate_pwelch(t, x, 'NumSegments', numSegs, 'MakeFigure', 0);
        else
            [f_tmp, sxx_tmp] = spectral_estimate_jlab(t, x, 'P', P, 'MakeFigure', 0);
        end

        %%% Skipping if segment wasn't long enough
        if isnan(f_tmp)
            sxx_all{i}(u,:) = NaN(size(f));
            continue
        end

        %%% Interpolating to common frequency grid for averaging
        idx = ~isnan(sxx_tmp);
        sxx_all{i}(u,:) = interp1(f_tmp(idx), sxx_tmp(idx), f, 'linear', NaN);

    end
    
end

%%% Plotting spectral estimate
hold on
f = f(:)';

for i = 1:length(moorings_all)
    mooring_no = moorings_all{i};
    x = sxx_all{i};

    mu = mean(log10(x), 1, 'omitnan');
    sd = std(log10(x), 0, 1, 'omitnan');

    lo = 10.^(mu - sd);
    hi = 10.^(mu + sd);

    % valid where both bounds are positive and finite (good for ylog)
    ok = isfinite(f) & isfinite(lo) & isfinite(hi) & (lo > 0) & (hi > 0);

    % Find contiguous valid segments
    d = diff([false, ok, false]);
    starts = find(d == 1);
    stops  = find(d == -1) - 1;

    % Draw one fill per segment (so it stops/restarts)
    for k = 1:numel(starts)
        idx = starts(k):stops(k);
        shaded_f = [f(idx), fliplr(f(idx))];
        shaded_x = [lo(idx), fliplr(hi(idx))];

        fill(shaded_f, shaded_x, colors(i,:), ...
            'FaceAlpha', 0.3, ...
            'EdgeColor', 'none', 'HandleVisibility', 'off');
    end

    plot(f, 10.^mu, 'Color', colors(i,:), 'LineWidth', lw, ...
        'DisplayName', string(upper(mooring_no)));
end

xlog; ylog;
xlim([1/180 3])
ylim([1e-7 1e0])
xlabel('Frequency [cycles/day]')
ylabel('Power Spectral Density')
grid on

text(0.5, 0.95, 'OBP', 'Units', 'Normalized', 'HorizontalAlignment', 'center', 'FontSize', fs, 'BackgroundColor', 'w');
text(0.5, 0.9, {'Unprocessed: Dashed', 'Detided: Solid'}, 'Units', 'Normalized', 'HorizontalAlignment', 'center', 'FontSize', fs-3)


%%% Adding legend
lgd = legend('Location','northeast');   % make legend first

% % Position legend in normalized figure units
% u = lgd.Units;
% lgd.Units = 'normalized';
% L = lgd.Position;          % [x y w h] in normalized figure coords
% lgd.Units = u;
% 
% % Add text box underneath the legend 
% annotation(gcf,'textbox', [L(1)-0.01, L(2)-0.1, L(3), 0.07], ...
%     'String', {'Unprocessed: Dashed', 'Detided: Solid'}, ...
%     'EdgeColor','none', ...
%     'FitBoxToText','on', ...
%     'HorizontalAlignment','center', ...
%     'VerticalAlignment','top');

%%% Adding tie points
xL = xlim;
yL = ylim;
y = 10^( log10(yL(1)) + 0.10*(log10(yL(2)) - log10(yL(1))) );

xline(1/100, 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Layer', 'bottom')
text(1/110, y, '100 days', 'Rotation', 90)
xline(1/20, 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Layer', 'bottom')
text(1/22.5, y, '20 days', 'Rotation', 90)
xline(1/3, 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Layer', 'bottom')
text(1/3.3, y, '3 days', 'Rotation', 90)
xline(1, 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Layer', 'bottom')
text(0.9, y, '1 day', 'Rotation', 90)
xline(2, 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Layer', 'bottom')
text(1.8, y, '12 hrs', 'Rotation', 90)

set(gca, 'FontSize', fs)

%%%%%%%%%%%%%%%%
%%% SLP data %%%
%%%%%%%%%%%%%%%%

nexttile()
hold on

%%% Computing spectral estimate for each deployment
clear sxx_all
for i = 1:length(moorings_all)
    mooring_no = moorings_all{i};

    for u = 1:length(bpr.(string(mooring_no) + '_start'))

        %%% Isolating annual segments
        t1 = bpr.(string(mooring_no) + '_start')(u);
        t2 = bpr.(string(mooring_no) + '_end')(u);
        idx = slp.datetime >= t1 & slp.datetime <= t2;
        t = slp.datetime(idx);
        x = slp.(string(mooring_no))(idx);

        %%% Skipping if segment is empty
        if isempty(x)
            sxx_all{i}(u,:) = NaN(size(f));
            continue
        end
    
        %%% Computing spectra
        if use_pwelch == 1
            [f_tmp, sxx_tmp] = spectral_estimate_pwelch(t, x,'NumSegments', numSegs, 'MakeFigure', 0);
        else
            [f_tmp, sxx_tmp] = spectral_estimate_jlab(t, x,'P', P, 'MakeFigure', 0);
        end

        %%% Skipping if segment wasn't long enough
        if isnan(f_tmp)
            sxx_all{i}(u,:) = NaN(size(f));
            continue
        end

        %%% Interpolating to common frequency grid for averaging
        idx = ~isnan(sxx_tmp);
        sxx_all{i}(u,:) = interp1(f_tmp(idx), sxx_tmp(idx), f);

    end
    
end

%%% Plotting spectral estimate
f = f(:)';

for i = 1:length(moorings_all)
    mooring_no = moorings_all{i};
    x = sxx_all{i};

    mu = mean(log10(x), 1, 'omitnan');
    sd = std(log10(x), 0, 1, 'omitnan');

    lo = 10.^(mu - sd);
    hi = 10.^(mu + sd);

    % valid where both bounds are positive and finite (good for ylog)
    ok = isfinite(f) & isfinite(lo) & isfinite(hi) & (lo > 0) & (hi > 0);

    % Find contiguous valid segments
    d = diff([false, ok, false]);
    starts = find(d == 1);
    stops  = find(d == -1) - 1;

    % Draw one fill per segment (so it stops/restarts)
    for k = 1:numel(starts)
        idx = starts(k):stops(k);
        shaded_f = [f(idx), fliplr(f(idx))];
        shaded_x = [lo(idx), fliplr(hi(idx))];

        fill(shaded_f, shaded_x, colors(i,:), ...
            'FaceAlpha', 0.3, ...
            'EdgeColor', 'none', 'HandleVisibility', 'off');
    end

    plot(f, 10.^mu, 'Color', colors(i,:), 'LineStyle', '-', 'LineWidth', lw, 'HandleVisibility', 'off');
end

xlog; ylog;
xlim([1/180 3])
ylim([1e-7 1e0])
yticklabels(' ')
xlabel('Frequency [cycles/day]')
grid on

text(0.5, 0.95, 'SLP', 'Units', 'Normalized', 'HorizontalAlignment', 'center', 'FontSize', fs, 'BackgroundColor', 'w');

%%% Adding tie points
xL = xlim;
yL = ylim;
y = 10^( log10(yL(1)) + 0.10*(log10(yL(2)) - log10(yL(1))) );

xline(1/100, 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Layer', 'bottom')
xline(1/20, 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Layer', 'bottom')
xline(1/3, 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Layer', 'bottom')
xline(1, 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Layer', 'bottom')
xline(2, 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Layer', 'bottom')
set(gca, 'FontSize', fs)

if use_pwelch == 1
    print(fig,'C:\Users\jak279\OneDrive - Yale University\Research\Figures\inverted_barometer_paper\bpr_slp_spectral_estimate_pwelch.png','-dpng', '-r600')
else 
    exportgraphics(fig,'C:\Users\jak279\OneDrive - Yale University\Research\Figures\inverted_barometer_paper\bpr_slp_spectral_estimate_jlab.png','Resolution',600)
end
