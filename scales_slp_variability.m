%% Script to show SLP coherence in the Arctic (SLP data kept in dbar)

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
    bpr.(string(mooring_no) + '_detided') = data.bpr.detided;
    bpr.(string(mooring_no) + '_start') = data.bpr.start;
    bpr.(string(mooring_no) + '_end') = data.bpr.end;

    %%% Loading SLP data
    data = load('C:\Users\jak279\OneDrive - Yale University\Research\Data\OBP_interp\SLP_mooring_' + string(mooring_no) + '_' + string(dt_string) + '.mat'); % SLP
    slp.datetime = data.slp.datetime;
   
    slp_tmp = data.slp.pressure; % basic slp, in dbar
    slp.(string(mooring_no)) = slp_tmp;

    % slp_tmp = -gsw_z_from_p(slp_tmp, data.slp.actual_lat);
    % slp.(string(mooring_no)) = NaN(size(data.slp.pressure));
    % for u = 1:length(bpr.(string(mooring_no) + '_start'))
    %     t1 = bpr.(string(mooring_no) + '_start')(u);
    %     t2 = bpr.(string(mooring_no) + '_end')(u);
    %     ind = data.slp.datetime >= t1 & data.slp.datetime <= t2;
    %     slp.(string(mooring_no))(ind) = slp_tmp(ind) - mean(slp_tmp(ind), 'omitnan'); % slp anomalies, in dbar
    % end
end

%% Loading gridded data

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

%% Applying Butterworth Filter

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

%%% Gridded SLP data
dt = hours(slp_grid.datetime(2) - slp_grid.datetime(1));
sampling_interval = minutes_per_hour * dt; % in minutes
slp_grid.highpassed = NaN(size(slp_grid.pressure));
slp_grid.lowpassed = NaN(size(slp_grid.pressure));
slp_grid.bandpassed = NaN(size(slp_grid.pressure));
for i = 1:size(slp_grid.pressure, 1)
    for j = 1:size(slp_grid.pressure, 2)
        [slp_grid.highpassed(i,j,:), ~] = highpassing_butterworth(slp_grid.pressure(i,j,:), sampling_interval, cutoff_lower);
        [tmp, slp_grid.lowpassed(i,j,:)] = highpassing_butterworth(slp_grid.pressure(i,j,:), sampling_interval, cutoff_upper);
        [~, slp_grid.bandpassed(i,j,:)] = highpassing_butterworth(tmp, sampling_interval, cutoff_lower);
    end
end

%% Plotting correlation between SLP at Mooring A and across the Arctic

mooring_no = 'a';
gridded = slp_grid;

%%% Beaufort Gyre region
lat1 = 80.5;
lat2 = 70.5;
lon1 = -170+360;
lon2 = -130+360;

fig = figure('Position', [10 10 1200 500]);
tiledlayout(1,3,'TileSpacing', 'none', 'Padding', 'none')
sgtitle('Relationship between SLP at Mooring ' + string(upper(mooring_no)) + ' and Arctic SLP')
fs = 12; lw = 2;

for m = 1:3
    if m == 1
        X = gridded.highpassed; % gridded time series
        Y = slp.(string(mooring_no) + '_highpassed'); % geographically isolated time series
        t = slp.datetime; % shared time grid
        title_str = '< 3 days';
    elseif m == 2
        X = gridded.bandpassed; % gridded time series
        Y = slp.(string(mooring_no) + '_bandpassed'); % geographically isolated time series
        t = slp.datetime; % shared time grid
        title_str = '3 - 20 days';
    elseif m == 3
        X = gridded.lowpassed; % gridded time series
        Y = slp.(string(mooring_no) + '_lowpassed'); % geographically isolated time series
        t = slp.datetime; % shared time grid
        title_str = '> 20 days';
    end

    %%% Conducting correlation analysis
    [Nx, Ny, Nt] = size(X);
    R = nan(Nx, Ny);
    p_val = nan(Nx, Ny);
    flagit = zeros(Nx, Ny);
    %%% Looping through each coordinate
    parfor i = 1:Nx
        for j = 1:Ny
            x = X(i,j,:);
            x = x(:);
            y = Y(:);
            y = y(:);

            %%% Finding good data
            good = isfinite(x) & isfinite(y);
            starts = find(diff([0; good]) == 1);
            if isempty(starts)
                continue
            end
            ends   = find(diff([good; 0]) == -1);

            [~, maxlength] = max(ends - starts);
            ind = starts(maxlength):ends(maxlength);
            x = x(ind);
            y = y(ind);

            %%% Computing significance levels
            if sum(isnan(x)) > 0.1*length(x)
                continue
            end
            rcrit = significance_testing_correlation_coefficient(x, y);

            %%% Calculating correlation coefficient
            [R(i,j), p_val(i,j)] = corr(x, y);

            %%% Flagging
            if abs(R(i,j)) < rcrit
                flagit(i,j) = 1;
            end
            % if p_val(i,j) > 0.05
            %     flagit(i,j) = 1;
            % end
        end
    end

    flagit = (flagit == 1);
    %R(flagit) = NaN;

    %%% Plotting
    nexttile()
    ax = axesm('stereo', 'Origin', [90 0], 'MapLatLimit', [50 90], 'Frame', 'on', 'Grid', 'on', 'MLineLocation', 60, 'MLabelLocation', 60, 'MLabelParallel', 48, 'MeridianLabel', 'on', 'PLineLocation', 15, 'ParallelLabel', 'on', 'MeridianLabel', 'on', 'LabelRotation', 'on');
    view(180, 90);
    pcolorm(gridded.lat, gridded.lon, R)
    land = readgeotable("landareas.shp");
    geoshow(ax, land, 'FaceColor', 'none')
    if m == 3
        cb = colorbar;
        ylabel(cb, {'Correlation Coefficient'}, 'FontSize', fs);
        cb.Label.Rotation = 270;
    end
    clim(1.*[-1 1])
    colormap(cmocean('balance'))

    % Overlay stippling where significant
    sig_mask = ~flagit;  % significant points
    sig_lat_pts = gridded.lat(sig_mask);
    sig_lon_pts = gridded.lon(sig_mask);
    scatterm(sig_lat_pts(:), sig_lon_pts(:), 1, 'k', 'filled', ...
        'MarkerFaceAlpha', 0.4)
    
    [~, c] = colornames('Crayola', 'white');
    textm([75, 78, 74, 90], [-150, -150, -140, 60], {'A', 'B', 'D', 'NP'}, 'FontSize', 8, 'FontWeight', 'b','Color', c)
    title(string(title_str))

    set(gca, 'FontSize', fs)
    set(ax, 'Color', 'none', 'Box', 'off', 'XColor', 'none', 'YColor', 'none');
end

exportgraphics(fig,'C:\Users\jak279\OneDrive - Yale University\Research\Figures\inverted_barometer_paper\scales_slp_variability.png','Resolution',600)

%% Plotting just synoptic scale correlation

mooring_no = 'a';
gridded = slp_grid;

%%% Creating figure
fig = figure('Position', [10 10 500 500]);
tiledlayout(1, 1, 'Padding', 'tight')
fs = 12; lw = 2;

%%% Settting data
X = gridded.bandpassed; % gridded time series
Y = slp.(string(mooring_no) + '_bandpassed'); % geographically isolated time series
t = slp.datetime; % shared time grid
title_str = 'Bandpassed 3-20 days';

%%% Conducting correlation analysis
[Nx, Ny, Nt] = size(X);
R = nan(Nx, Ny);
p_val = nan(Nx, Ny);
flagit = zeros(Nx, Ny);
%%% Looping through each coordinate
for i = 1:Nx
    for j = 1:Ny
        x = X(i,j,:);
        x = x(:);
        y = Y(:);
        y = y(:);

        %%% Finding good data
        good = isfinite(x) & isfinite(y);
        starts = find(diff([0; good]) == 1);
        if isempty(starts)
            continue
        end
        ends   = find(diff([good; 0]) == -1);

        [~, maxlength] = max(ends - starts);
        ind = starts(maxlength):ends(maxlength);
        x = x(ind);
        y = y(ind);

        %%% Computing significance levels
        if sum(isnan(x)) > 0.1*length(x)
            continue
        end
        rcrit = significance_testing_correlation_coefficient(x, y);

        %%% Calculating correlation coefficient
        [R(i,j), p_val(i,j)] = corr(x, y);

        %%% Flagging
        if abs(R(i,j)) < rcrit
            flagit(i,j) = 1;
        end
        % if p_val(i,j) > 0.05
        %     flagit(i,j) = 1;
        % end
    end
end

flagit = (flagit == 1);
%R(flagit) = NaN;

%%% Plotting
ax = axesm('stereo', 'Origin', [90 0], 'MapLatLimit', [59.99 90], 'Frame', 'on', 'Grid', 'on', 'MLineLocation', 60, 'MLabelLocation', 60, 'MLabelParallel', 'south', 'MeridianLabel', 'on', 'PLineLocation', [65 75 85], 'ParallelLabel', 'on', 'MeridianLabel', 'on', 'LabelRotation', 'on');
view(180, 90);

lon_wrap  = [gridded.lon;  gridded.lon(1,:) + 360];
lat_wrap  = [gridded.lat;  gridded.lat(1,:)];
data_mean = R;
data_wrap = [data_mean;     data_mean(1,:)];
surfm(lat_wrap', lon_wrap', data_wrap');

land = readgeotable("landareas.shp");
geoshow(ax, land, 'FaceColor', 'none')
cb = colorbar;
ylabel(cb, {'Correlation Coefficient'}, 'FontSize', fs);
cb.Label.Rotation = 270;
clim(1.*[-1 1])
colormap(cmocean('balance'))

% Overlay stippling where significant
sig_mask = ~flagit;  % significant points
sig_lat_pts = gridded.lat(sig_mask);
sig_lon_pts = gridded.lon(sig_mask);
scatterm(sig_lat_pts(:), sig_lon_pts(:), 1, 'k', 'filled', ...
    'MarkerFaceAlpha', 0.4)

[~, c] = colornames('Crayola', 'white');
textm([75, 78, 74, 90], [-150, -150, -140, 60], {'A', 'B', 'D', 'NP'}, 'FontSize', fs, 'FontWeight', 'b','Color', c, 'HorizontalAlignment', 'center')
%title(string(title_str))

set(gca, 'FontSize', fs)
set(ax, 'Color', 'none', 'Box', 'off', 'XColor', 'none', 'YColor', 'none');

exportgraphics(fig,'C:\Users\jak279\OneDrive - Yale University\Research\Figures\inverted_barometer_paper\scales_slp_variability_synoptic.png','Resolution',600)

