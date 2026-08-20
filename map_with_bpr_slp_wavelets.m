%% Script to visualize bpr/slp data
clc; clear;

%% Loading mooring data

moorings_all = {'a', 'b', 'd', 'np'};
dt_string = '6hr';

for i = 1:length(moorings_all)

    %%% Selecting mooring
    mooring_no = moorings_all{i};

    %%% Loading BPR data
    data = load('C:\Users\jak279\OneDrive - Yale University\Research\Data\OBP_interp\BPR_mooring_' + string(mooring_no) + '_' + string(dt_string) + '.mat'); % BPR
    bpr.datetime = data.bpr.datetime;
    bpr_tmp = data.bpr.detided;
    bpr.(string(mooring_no) + '_detided') = bpr_tmp;
    bpr.(string(mooring_no) + '_start') = data.bpr.start;
    bpr.(string(mooring_no) + '_end') = data.bpr.end;

    %%% Loading SLP data
    data = load('C:\Users\jak279\OneDrive - Yale University\Research\Data\OBP_interp\SLP_mooring_' + string(mooring_no) + '_' + string(dt_string) + '.mat'); % SLP
    slp.datetime = data.slp.datetime;
    slp_tmp = data.slp.pressure;
    slp.(string(mooring_no)) = slp_tmp;
    %slp_tmp = -gsw_z_from_p(slp_tmp, data.slp.actual_lat);
    slp.(string(mooring_no)) = NaN(size(slp_tmp));
    for u = 1:length(bpr.(string(mooring_no) + '_start'))
        t1 = bpr.(string(mooring_no) + '_start')(u);
        t2 = bpr.(string(mooring_no) + '_end')(u);
        ind = data.slp.datetime >= t1 & data.slp.datetime <= t2;
        slp.(string(mooring_no))(ind) = slp_tmp(ind) - mean(slp_tmp(ind), 'omitnan');
    end
end


%% Loading gridded slp data

%%% Years of interest
yrs = 2003:2025;

%%% SLP
slp_full = loadingGriddedSLPdata(yrs);

%%% Getting break points for anomaly calculation
mean_start = mean([bpr.a_start; bpr.b_start; bpr.d_start], 1, 'omitnan');
mean_end = mean([bpr.a_end; bpr.b_end; bpr.d_end], 1, 'omitnan');

%%% Computing anomalies
slp_full.pressure_anom = NaN(size(slp_full.pressure));
for i = 1:length(mean_start)
    ind = slp_full.datetime >= mean_start(i) & slp_full.datetime <= mean_end(i);
    slp_full.pressure_anom(:,:,ind) = slp_full.pressure(:,:,ind) - mean(slp_full.pressure(:,:,ind), 3, 'omitnan');
end

%%% Running mean to get slp on chosen time grid
clear slp_grid
[slp_grid.datetime, slp_grid.pressure] =  running_mean(slp_full.datetime, slp_full.pressure, bpr.datetime);
[slp_grid.datetime, slp_grid.pressure_anom] = running_mean(slp_full.datetime, slp_full.pressure_anom, bpr.datetime);
slp_grid.lat = slp_full.lat;
slp_grid.lon = slp_full.lon;
clear slp_full

%% Bandpassing data

minutes_per_hour = 60;
hours_per_day = 24;
cutoff_lower = hours_per_day * 3; % in hours
cutoff_upper = hours_per_day * 20; % in hours

%%% OBP at mooring sites
dt = hours(bpr.datetime(2) - bpr.datetime(1));
sampling_interval = minutes_per_hour * dt; % in minutes
for i = 1:length(moorings_all)
    mooring_no = moorings_all(i);
    [bpr.(string(mooring_no) + '_highpassed'), ~] = highpassing_butterworth(bpr.(string(mooring_no) + '_detided'), sampling_interval, cutoff_lower);
    [tmp, bpr.(string(mooring_no) + '_lowpassed')] = highpassing_butterworth(bpr.(string(mooring_no) + '_detided'), sampling_interval, cutoff_upper);
    [~, bpr.(string(mooring_no) + '_bandpassed')] = highpassing_butterworth(tmp, sampling_interval, cutoff_lower);
end

%%% SLP at mooring sites
dt = hours(slp.datetime(2) - slp.datetime(1));
sampling_interval = minutes_per_hour * dt; % in minutes
for i = 1:length(moorings_all)
    mooring_no = moorings_all(i);
    [slp.(string(mooring_no) + '_highpassed'), ~] = highpassing_butterworth(slp.(string(mooring_no)), sampling_interval, cutoff_lower);
    [tmp, slp.(string(mooring_no) + '_lowpassed')] = highpassing_butterworth(slp.(string(mooring_no)), sampling_interval, cutoff_upper);
    [~, slp.(string(mooring_no) + '_bandpassed')] = highpassing_butterworth(tmp, sampling_interval, cutoff_lower);
end


%%

% %%% Loading bathymetry data
% if ~exist('bathylats', 'var')
%     %%% Loading bathymetry data
%     [A,R] = readgeoraster("IBCAO_V3_500m_SM.tif","OutputType","double");
%     % Create grid of X,Y values
%     [x,y] = worldGrid(R);
%     % Convert grid of X,Y values to latitude/longitude
%     [bathylats,bathylons] = projinv(R.ProjectedCRS,x,y);
%     A(bathylats <= 60) = NaN;
%     A = -A;
% end

%%% Making figure
fig = figure('Position', [10 10 1400 600]);
tiledlayout(3,2, 'Padding', 'tight', 'TileSpacing', 'tight')
fs = 12; lw = 2;

%%% Setting mooring
mooring_no = 'a';
idx = 20;
t1 = bpr.(string(mooring_no) + '_start')(idx);
t2 = bpr.(string(mooring_no) + '_end')(idx);

%%%%%%%%%%%%%%%%%%%
%%% Map subplot %%%
%%%%%%%%%%%%%%%%%%%
nexttile(1, [3, 1])
ax = axesm('stereo', 'Origin', [90 0], 'MapLatLimit', [59.99 90], 'Frame', 'on', 'Grid', 'on', 'MLineLocation', 60, 'MLabelLocation', 60, 'MLabelParallel', 'south', 'MeridianLabel', 'on', 'PLineLocation', [65 75 85], 'ParallelLabel', 'on', 'MeridianLabel', 'on', 'LabelRotation', 'on');
view(180, 90);

%%% Plotting bathymetry
% delta_bathy = 4;
% pcolorm(bathylats(1:delta_bathy:11617, 1:delta_bathy:11617), bathylons(1:delta_bathy:11617, 1:delta_bathy:11617), A(1:delta_bathy:11617, 1:delta_bathy:11617));
% colormap(cmocean('deep')); clim([0 5000]);

%%% Plotting SLP
ind = slp_grid.datetime >= datetime(2025, 1, 1, 'TimeZone', 'UTC') & slp_grid.datetime < datetime(2026, 1, 1, 'TimeZone', 'UTC');
lon_wrap  = [slp_grid.lon;  slp_grid.lon(1,:) + 360];
lat_wrap  = [slp_grid.lat;  slp_grid.lat(1,:)];
data_mean = mean(slp_grid.pressure(:,:,ind), 3, 'omitnan');
data_wrap = [data_mean;     data_mean(1,:)];
surfm(lat_wrap', lon_wrap', data_wrap');
shading flat
colormap(cmocean('balance')); clim([10.02 10.18]); 
cb = colorbar;
cb.Ticks = 10:0.04:10.2;
title('Sea Level Pressure [dbar]', 'FontSize', fs, 'FontWeight', 'normal')

land = readgeotable("landareas.shp");
geoshow(ax, land, 'FaceColor', [0.8 0.8 0.8])
textm([75, 78, 74, 90], [-150, -150, -140, 60], {'A', 'B', 'D', 'NP'}, 'FontSize', 12, 'FontWeight', 'b','Color', 'k', 'HorizontalAlignment', 'center')
set(gca, 'FontSize', fs)
set(ax, 'Color', 'none') 
box off
set(ax, 'XColor', 'none')
set(ax, 'YColor', 'none')
text(0.2, 0.9, '(a)', 'Units', 'normalized', 'Color', 'r',...
    'HorizontalAlignment', 'center', 'FontSize', fs, ...
    'FontWeight', 'bold', 'BackgroundColor', 'none', 'Parent', ax)

%%%%%%%%%%%%%%%%%%%
%%% OBP Wavelet %%%
%%%%%%%%%%%%%%%%%%%

c1 = [1 0.6 0];

ax(1) = nexttile(2);
x = bpr.(string(mooring_no) + '_detided');
t = bpr.datetime;
dt = datenum(t(2)) - datenum(t(1));
fs_num = 1./dt;
plot_wavelet(ax(1), t, x, fs_num, 'CLim', log10([0 4]),'YTicks', [1 4 16 64 256 1024]);
xline(yearfrac(datenum([t1 t2])), 'Color', c1, 'LineWidth', lw)
xticklabels(' ')
text(0.5, 0.93, 'OBP', 'Units', 'normalized', 'Color', 'r',...
    'HorizontalAlignment', 'center', 'FontSize', fs, ...
    'FontWeight', 'bold', 'BackgroundColor', 'none', 'Parent', ax(1))
text(0.03, 0.93, '(b)', 'Units', 'normalized', 'Color', 'r',...
    'HorizontalAlignment', 'center', 'FontSize', fs, ...
    'FontWeight', 'bold', 'BackgroundColor', 'none', 'Parent', ax(1))

%%%%%%%%%%%%%%%%%%%
%%% SLP Wavelet %%%
%%%%%%%%%%%%%%%%%%%

ax(2) = nexttile(4);
x = slp.(string(mooring_no));
t = slp.datetime;
dt = datenum(t(2)) - datenum(t(1));
fs_num = 1./dt;
plot_wavelet(ax(2), t, x, fs_num, 'Clim', log10([0 4]), 'YTicks', [1 4 16 64 256 1024]);
xline(yearfrac(datenum([t1 t2])), 'Color', c1, 'LineWidth', lw)
text(0.5, 0.93, 'SLP', 'Units', 'normalized', 'Color', 'r',...
    'HorizontalAlignment', 'center', 'FontSize', fs, ...
    'FontWeight', 'bold', 'BackgroundColor', 'none', 'Parent', ax(2))
text(0.03, 0.93, '(c)', 'Units', 'normalized', 'Color', 'r',...
    'HorizontalAlignment', 'center', 'FontSize', fs, ...
    'FontWeight', 'bold', 'BackgroundColor', 'none', 'Parent', ax(2))

%%%%%%%%%%%%%%%%%%%
%%% Time Series %%%
%%%%%%%%%%%%%%%%%%%

ax(3) = nexttile(6);
yyaxis left
hold on
plot(bpr.datetime, bpr.(string(mooring_no) + '_detided'), 'LineWidth', lw, 'DisplayName', 'OBP')
ylabel({'OBP Anomaly', '[dbar]'})
ylim(0.2.*[-1 1])
yticks(-0.1:0.1:0.1)
yyaxis right
hold on
plot(slp.datetime, slp.(string(mooring_no)), 'LineWidth', lw, 'DisplayName', 'SLP')
ylabel({'SLP Anomaly', '[dbar]'}, 'Rotation', 270)
ylim(0.4.*[-1 1])
yticks(-0.2:0.2:0.2)
xlim([t1 t2])
grid on
yline(0, '-', 'Color', [0.5 0.5 0.5], 'LineWidth', lw)
text(0.03, 0.93, '(d)', 'Units', 'normalized', 'Color', 'r',...
    'HorizontalAlignment', 'center', 'FontSize', fs, ...
    'FontWeight', 'bold', 'BackgroundColor', 'none', 'Parent', ax(3))

set(gca, 'FontSize', fs)

% Get positions of the axes in figure coordinates [x y width height]
pos_slp  = ax(2).Position;   % SLP wavelet
pos_ts   = ax(3).Position;   % Time series

% Get the x-limits of the SLP wavelet to convert t1,t2 to normalized coords
xlim_slp = ax(2).XLim;
t1_norm  = yearfrac(datenum(t1));
t2_norm  = yearfrac(datenum(t2));

% Convert t1 and t2 to figure x coordinates within the SLP wavelet axes
x1_fig = pos_slp(1) + (t1_norm - xlim_slp(1)) / (xlim_slp(2) - xlim_slp(1)) * pos_slp(3);
x2_fig = pos_slp(1) + (t2_norm - xlim_slp(1)) / (xlim_slp(2) - xlim_slp(1)) * pos_slp(3);

% Bottom of SLP wavelet axes
y_slp_bottom = pos_slp(2);

% Top of time series axes
y_ts_top = pos_ts(2) + pos_ts(4);

% Left edge of time series axes (t1)
x1_ts_fig = pos_ts(1);
% Right edge of time series axes (t2)
x2_ts_fig = pos_ts(1) + pos_ts(3);

% Create invisible overlay axes covering the whole figure
ax_overlay = axes('Position', [0 0 1 1], 'Visible', 'off');
ax_overlay.XLim = [0 1];
ax_overlay.YLim = [0 1];
hold(ax_overlay, 'on')

% Draw the connecting lines with alpha
l1 = plot(ax_overlay, [x1_fig, x1_ts_fig], [y_slp_bottom, y_ts_top], ...
    '-', 'Color', c1, 'LineWidth', lw);
l2 = plot(ax_overlay, [x2_fig, x2_ts_fig], [y_slp_bottom, y_ts_top], ...
    '-', 'Color', c1, 'LineWidth', lw);

% Set alpha
l1.Color(4) = 0.4;   % 0 = fully transparent, 1 = fully opaque
l2.Color(4) = 0.4;

% Make sure overlay is on top but doesn't interfere
uistack(ax_overlay, 'top')

%%% Saving figure
print(fig,'C:\Users\jak279\OneDrive - Yale University\Research\Figures\inverted_barometer_paper\map_bpr_slp_wavelets', '-dpng', '-r600')

%% Figure showing cross-power and bandpassed

fig = figure('Position', [10 10 1200 600]);
tiledlayout(2, 1, 'Padding', 'tight', 'TileSpacing', 'tight')
clear ax

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Wavelet Cross-Power %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

ind = slp.datetime >= bpr.datetime(1) & slp.datetime <= bpr.datetime(end);
tx = slp.datetime(ind);
x = slp.(string(mooring_no))(ind);
ty = bpr.datetime;
y = bpr.(string(mooring_no) + '_detided');

%%% Wavelet coherence
ax(1) = nexttile();
plot_wavelet_coherence(ax(1), [], tx, x, y, 'PlotPhase', false, 'RunMonteCarlo', true, ...
    'NSurrogates', 200, 'SigLevel', 95, 'YTicks', [1 4 16 64 256 1024]);
hold on
xline(yearfrac(datenum([t1 t2])), 'Color', c1, 'LineWidth', lw)
xlabel(' ')
set(gca, 'FontSize', fs)
text(0.02, 0.93, '(a)', 'Units', 'normalized', 'Color', 'r',...
    'HorizontalAlignment', 'center', 'FontSize', fs, ...
    'FontWeight', 'bold', 'BackgroundColor', 'none', 'Parent', ax(1))


%%%%%%%%%%%%%%%%%%%%%%%
%%% Bandpassed Data %%%
%%%%%%%%%%%%%%%%%%%%%%%

ax(2) = nexttile();
yyaxis left
hold on
plot(bpr.datetime, bpr.(string(mooring_no) + '_bandpassed'), 'LineWidth', lw)
ylabel({'OBP Anomaly', '[dbar]'})
ylim(0.15.*[-1 1])
yyaxis right
hold on
plot(slp.datetime, slp.(string(mooring_no) + '_bandpassed'), 'LineWidth', lw)
ylabel({'SLP Anomaly', '[dbar]'}, 'Rotation', 270)
ylim(0.3.*[-1 1])
grid on
xlim([t1 t2])
set(gca, 'FontSize', fs)
yline(0, 'Color', [0.5 0.5 0.5], 'LineWidth', lw)
text(0.5, 0.93, 'Bandpassed 3-20 days', 'Units', 'normalized', ...
    'HorizontalAlignment', 'center', 'FontSize', fs, ...
    'FontWeight', 'bold', 'BackgroundColor', 'none', 'Parent', ax(2))
text(0.02, 0.93, '(b)', 'Units', 'normalized', 'Color', 'r',...
    'HorizontalAlignment', 'center', 'FontSize', fs, ...
    'FontWeight', 'bold', 'BackgroundColor', 'none', 'Parent', ax(2))

% Get positions of the axes in figure coordinates [x y width height]
pos_slp  = ax(1).Position;   % BPR-SLP cross-wavelet
pos_ts   = ax(2).Position;   % Time series

% Get the x-limits of the BPR-SLP cross-wavelet to convert t1,t2 to normalized coords
xlim_slp = ax(1).XLim;
t1_norm  = yearfrac(datenum(t1));
t2_norm  = yearfrac(datenum(t2));

% Convert t1 and t2 to figure x coordinates within the SLP wavelet axes
x1_fig = pos_slp(1) + (t1_norm - xlim_slp(1)) / (xlim_slp(2) - xlim_slp(1)) * pos_slp(3);
x2_fig = pos_slp(1) + (t2_norm - xlim_slp(1)) / (xlim_slp(2) - xlim_slp(1)) * pos_slp(3);

% Bottom of SLP wavelet axes
y_slp_bottom = pos_slp(2);

% Top of time series axes
y_ts_top = pos_ts(2) + pos_ts(4);

% Left edge of time series axes (t1)
x1_ts_fig = pos_ts(1);
% Right edge of time series axes (t2)
x2_ts_fig = pos_ts(1) + pos_ts(3);

% Create invisible overlay axes covering the whole figure
ax_overlay = axes('Position', [0 0 1 1], 'Visible', 'off');
ax_overlay.XLim = [0 1];
ax_overlay.YLim = [0 1];
hold(ax_overlay, 'on')

% Draw the connecting lines with alpha
l1 = plot(ax_overlay, [x1_fig, x1_ts_fig], [y_slp_bottom, y_ts_top], ...
    '-', 'Color', c1, 'LineWidth', lw);
l2 = plot(ax_overlay, [x2_fig, x2_ts_fig], [y_slp_bottom, y_ts_top], ...
    '-', 'Color', c1, 'LineWidth', lw);

% Set alpha
l1.Color(4) = 0.4;   % 0 = fully transparent, 1 = fully opaque
l2.Color(4) = 0.4;


%%% Saving figure
print(fig,'C:\Users\jak279\OneDrive - Yale University\Research\Figures\inverted_barometer_paper\bpr_slp_crosspower_bandpassed', '-dpng', '-r600')