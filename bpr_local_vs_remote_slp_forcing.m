%% Script to visualize coherence/phase between OBP and local and remote SLP

clc; clear;

%% Loading mooring data

moorings_all = {'a', 'b', 'd', 'np'};
dt_string = '1day';

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
    bpr.(string(mooring_no) + '_lat') = data.bpr.lat;
    lon = data.bpr.lon;
    if lon < 0
        lon = lon + 360;
    end
    bpr.(string(mooring_no) + '_lon') = lon;

    %%% Loading SLP data
    data = load('C:\Users\jak279\OneDrive - Yale University\Research\Data\OBP_interp\SLP_mooring_' + string(mooring_no) + '_' + string(dt_string) + '.mat'); % SLP
    slp.datetime = data.slp.datetime;
    slp_tmp = data.slp.pressure;
    slp.(string(mooring_no)) = slp_tmp;
    % slp.(string(mooring_no)) = NaN(size(slp_tmp));
    % for u = 1:length(bpr.(string(mooring_no) + '_start'))
    %     t1 = bpr.(string(mooring_no) + '_start')(u);
    %     t2 = bpr.(string(mooring_no) + '_end')(u);
    %     ind = data.slp.datetime >= t1 & data.slp.datetime <= t2;
    %     slp.(string(mooring_no))(ind) = slp_tmp(ind) - mean(slp_tmp(ind), 'omitnan');
    % end
end

clear slp_tmp ind t1 t2 u data bpr_tmp i 

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

%% Loading gridded wind data

%%% Years of interest
yrs = 2003:2025;

%%% Wind
wind_full = loadingGriddedWinds(yrs, 'NCEP');

%%% Running mean to get wind on chosen time grid
clear wind_grid
[wind_grid.datetime, wind_grid.u_speed] =  running_mean(wind_full.datetime, wind_full.u_speed, bpr.datetime);
[wind_grid.datetime, wind_grid.v_speed] = running_mean(wind_full.datetime, wind_full.v_speed, bpr.datetime);
[wind_grid.datetime, wind_grid.curl] = running_mean(wind_full.datetime, wind_full.curl, bpr.datetime);
wind_grid.lat = wind_full.lat;
wind_grid.lon = wind_full.lon;
clear wind_full

%% Computing coherence/phase for key coordinates

mooring_no = 'd';
key_latitudes   = [74, 78, 82, 86, 90, 87, 84, 81, 78, 75, 72, bpr.(string(mooring_no) + '_lat')];
key_longitudes   = [215, 215, 215, 215, 0, 10, 15, 10, 5, 0, 355, bpr.(string(mooring_no) + '_lon')];
names = {'Canada Basin', 'Canada Basin N1', 'Canada Basin N2', 'Canada Basin N3', 'North Pole', 'Nansen Basin N', 'Nansen Basin', 'Nansen Basin S', 'Yermak Plateau', 'Fram Strait', 'Norwegian Sea', 'Mooring ' + string(upper(mooring_no))};

%%% Pre-allocating cells to hold the spectral analysis
[coherence_all, cospectrum_all, transfer_function_all, ...
    phase_all, theor_coherence_crit_all, coherence_crit_all,...
    cospectrum_crit_all, transfer_function_crit_all] = deal(cell(1, length(key_longitudes)));

numSegs = 8;
numSurr = 50;

%%% Looping through key coordinates
for i = 1:length(key_latitudes)

    lat = key_latitudes(i);
    lon = key_longitudes(i);
    label = names{i};

    %%% Extracting slp data at key coordinate
    isLat = slp_grid.lat > lat-3 & slp_grid.lat < lat+3; % find latitude indices around target point
    isLon = slp_grid.lon > lon-3 & slp_grid.lon < lon+3; % find longitude indices around target point
    nearI = find(isLat & isLon);
    nNear = length(nearI); % all the indices in the vicinity
    dist = distance(slp_grid.lat(nearI),slp_grid.lon(nearI),repmat(lat,nNear,1),repmat(lon,nNear,1)); % compute distance from all the indices in the vicinity
    [~,I] = min(dist);
    closestI = nearI(I);
    [r,c] = ind2sub(size(slp_grid.lat), closestI); % index of the closest grid cell

    %%% Creating slp time series 
    slp_key_coordinate = squeeze(slp_grid.pressure(r,c,:));

    %%% Cross-spectral analysis
    [f, coherence_all{i}, cospectrum_all{i}, transfer_function_all{i}, phase_all{i}, theor_coherence_crit_all{i}, coherence_crit_all{i}, cospectrum_crit_all{i}, transfer_function_crit_all{i}] = coherence_phase_pwelch_full_record(bpr.datetime, bpr.(string(mooring_no) + '_detided'), slp_grid.datetime, slp_key_coordinate, 'StartDates', bpr.(string(mooring_no) + '_start'), 'EndDates', bpr.(string(mooring_no) + '_end'), 'NumSegments', numSegs, 'Nsurr', numSurr, 'MonteCarloSignificance', 1);

end

%%% Plotting
fig = figure('Position', [10 10 1200 400]);
tiledlayout(1, 4, 'TileSpacing', 'tight', 'Padding', 'none');
colors = cmocean('thermal', length(key_longitudes)-1);
colors = vertcat(colors, [0,1,0]);
%colors = orderedcolors('gem12');
sgtitle('Relationship with SLP at Mooring ' + string(upper(mooring_no)));
lw = 2;

%%% Coherence
ax = nexttile();
plot_mean_sig_thresh(ax, f, coherence_all, coherence_crit_all, names, colors, lw, ...
    'YLim',[0 1], 'YTicks',0:0.25:1, 'XLog',true, 'YLog',false, ...
    'YLabel', 'Magnitude-Squared Coherence', ...
    'XLabel', 'Frequency [cycles/day]',...
    'ShowLegend', 0, 'LegendLocation','northwest', ...
    'ShowXTicks', true, ...
    'TiePeriods',[100 20 3], 'TieLabels',{'100 days','20 days','3 days'}, 'ShowTieText', true);

%%% Cross-Spectrum
ax = nexttile();
plot_mean_sig_thresh(ax, f, cospectrum_all, cospectrum_crit_all, names, colors, lw, ...
    'YLim',[0 0.5], 'XLog',true, 'YLog',true, ...
    'YLabel', "Cross-Spectrum", ...
    'XLabel', 'Frequency [cycles/day]',...
    'ShowLegend', 0, 'LegendLocation', 'southwest', ...
    'ShowXTicks', true, 'TiePeriods',[100 20 3], 'ShowTieText', false);

%%% Transfer Function
ax = nexttile();
plot_mean_sig_thresh(ax, f, transfer_function_all, transfer_function_crit_all, names, colors, lw, ...
    'YLim',1.*[0 1], 'YTicks',-1:0.25:1, 'XLog',true, 'YLog',false, ...
    'YLabel', "Transfer Function Magnitude", ...
    'XLabel', 'Frequency [cycles/day]',...
    'ShowXTicks', true, 'TiePeriods',[100 20 3], 'ShowTieText', false, 'isTransferFunction', 0);
yline(ax, 0, '-', 'Color', [0.5 0.5 0.5], 'LineWidth', lw);

%%% Phase
ax = nexttile();
phase_all_corrected = phase_all;
for i = 1:size(phase_all_corrected, 2)
    ind = coherence_all{i} < coherence_crit_all{i};
    phase_all_corrected{i}(ind) = NaN;
    phase_all_corrected{i} = rad2deg(unwrap(phase_all_corrected{i}));
end
plot_mean_sig_thresh(ax, f, phase_all_corrected, cell(size(phase_all_corrected)), names, colors, lw, ...
    'YLim',[-180 180], 'YTicks',-180:60:180, 'XLog',true, 'YLog',false, ...
    'XLabel',"Frequency [cycles/day]", ...
    'YLabel', "Phase [deg]", ...
    'ShowLegend', 1, 'LegendLocation', 'eastoutside',...
    'ShowXTicks', true, 'TiePeriods',[100 20 3], 'ShowTieText', false);
yline(ax, 0, '-', 'Color', [0.5 0.5 0.5], 'LineWidth', lw, 'HandleVisibility', 'off');

%% Map showing "key coordinates"

%%% Loading bathymetry data
if ~exist('bathylats', 'var')
    %%% Loading bathymetry data
    [A,R] = readgeoraster("IBCAO_V3_500m_SM.tif","OutputType","double");
    % Create grid of X,Y values
    [x,y] = worldGrid(R);
    % Convert grid of X,Y values to latitude/longitude
    [bathylats,bathylons] = projinv(R.ProjectedCRS,x,y);
    A(bathylats < 60) = NaN;
    A = -A;
end

%%% Making figure
figure()
fs = 12;
ax = axesm('stereo', 'Origin', [90 0], 'MapLatLimit', [60 90], 'Frame', 'on', 'Grid', 'on', 'MLineLocation', 60, 'MLabelLocation', 60, 'MLabelParallel', 48, 'MeridianLabel', 'on', 'PLineLocation', 15, 'ParallelLabel', 'on', 'MeridianLabel', 'on', 'LabelRotation', 'on');
view(180, 90);
pcolorm(bathylats(1:4:11617, 1:4:11617), bathylons(1:4:11617, 1:4:11617), A(1:4:11617, 1:4:11617));
colormap(cmocean('deep')); clim([0 5000]);
scatterm(key_latitudes, key_longitudes, 50, colors, 'filled')
land = readgeotable("landareas.shp");
geoshow(ax, land, 'FaceColor', [0.8 0.8 0.8])

textm([75, 78, 74, 90], [-150, -150, -140, 60], {'A', 'B', 'D', 'NP'}, 'FontSize', 8, 'FontWeight', 'b','Color', 'k')

set(gca, 'FontSize', fs)
set(ax, 'Color', 'none', 'Box', 'off', 'XColor', 'none', 'YColor', 'none');

%% Looping through SLP grid and computing average coherence/phase in desired frequency band

mooring_no = 'np'; 
gridded = wind_grid;

%%% Pre-allocating cells to hold the spectral analysis
[coherence_all, cospectrum_all, transfer_function_all, ...
    phase_all, theor_coherence_crit_all, coherence_crit_all,...
    cospectrum_crit_all, transfer_function_crit_all] = deal(cell(size(gridded.lat)));
coherence_map = NaN(size(gridded.lat));
phase_map = NaN(size(gridded.lat));

numSegs = 8;
numSurr = 50;

for i = 1:1:size(gridded.lat, 1)
    disp(i)
    for j = 1:1:size(gridded.lat, 2)
        if gridded.lat(i,j) < 60 
            continue
        end
        %%% Creating time series
        ind = gridded.datetime < datetime(2010, 1, 1, 'TimeZone', 'UTC');
        gridded_time = gridded.datetime(ind);
        gridded_coordinate = squeeze(gridded.v_speed(i,j,ind));

        %%% Geographically isolated time series to compare to
        ind = bpr.datetime < datetime(2010, 1, 1, 'TimeZone', 'UTC');
        t = bpr.datetime(ind);
        x = bpr.(string(mooring_no) + '_detided');
        x = x (ind);

        %%% Cross-spectral analysis
        [f, coherence_all{i,j}, cospectrum_all{i,j}, transfer_function_all{i,j}, phase_all{i,j}, theor_coherence_crit_all{i,j}, coherence_crit_all{i,j}, cospectrum_crit_all{i,j}, transfer_function_crit_all{i,j}] = coherence_phase_pwelch_full_record(t, x, gridded_time, gridded_coordinate, 'StartDates', bpr.(string(mooring_no) + '_start'), 'EndDates', bpr.(string(mooring_no) + '_end'), 'NumSegments', numSegs, 'Nsurr', numSurr, 'MonteCarloSignificance', 0);
        ind = f >= 1/30;
        coherence_map(i,j) = mean(coherence_all{i,j}(ind), 'omitnan');
        phase_all_corrected = phase_all{i,j};
        idx = coherence_all{i,j} < 0.1;
        phase_all_corrected(idx) = NaN;
        phase_all_corrected = rad2deg(unwrap(phase_all_corrected));
        phase_map(i,j) = mean(phase_all_corrected(ind), 'omitnan');
    end
end

%% Plotting map of coherence/phase in desired frequency band

figure('Position', [10 10 1000 500])
tiledlayout(1,2, 'TileSpacing', 'tight')
sgtitle('Relationship between OBP at Mooring ' + string(upper(mooring_no)) + ' and Arctic meridional winds averaged over the < 30 days frequency band')

nexttile()
fs = 12;
gridded = wind_grid;
ax1 = axesm('stereo', 'Origin', [90 0], 'MapLatLimit', [59.99 90], 'Frame', 'on', 'Grid', 'on', 'MLineLocation', 60, 'MLabelLocation', 60, 'MLabelParallel', 60, 'MeridianLabel', 'on', 'PLineLocation', [65 75 85], 'ParallelLabel', 'on', 'MeridianLabel', 'on', 'LabelRotation', 'on');
view(180, 90);
pcolorm(gridded.lat, gridded.lon, coherence_map)
land = readgeotable("landareas.shp");
geoshow(ax1, land, 'FaceColor', 'none')
clim(0.5.*[0 1])
cb = colorbar;
ylabel(cb, 'Magnitude-Squared Coherence', 'Rotation', 270)
colormap(ax1, 'hot')
[~, c] = colornames('Crayola', 'white');
textm([75, 78, 74, 90], [-150, -150, -140, 60], {'A', 'B', 'D', 'NP'}, 'FontSize', fs, 'FontWeight', 'b','Color', c, 'HorizontalAlignment', 'center')
set(gca, 'FontSize', fs)
set(ax1, 'Color', 'none', 'Box', 'off', 'XColor', 'none', 'YColor', 'none');

nexttile()
ax2 = axesm('stereo', 'Origin', [90 0], 'MapLatLimit', [59.99 90], 'Frame', 'on', 'Grid', 'on', 'MLineLocation', 60, 'MLabelLocation', 60, 'MLabelParallel', 60, 'MeridianLabel', 'on', 'PLineLocation', [65 75 85], 'ParallelLabel', 'on', 'MeridianLabel', 'on', 'LabelRotation', 'on');
view(180, 90);
pcolorm(gridded.lat, gridded.lon, phase_map)
land = readgeotable("landareas.shp");
geoshow(ax2, land, 'FaceColor', 'none')
clim([-180 180])
cb = colorbar;
ylabel(cb, 'Phase', 'Rotation', 270)
colormap(ax2, 'hsv')
[~, c] = colornames('Crayola', 'white');
textm([75, 78, 74, 90], [-150, -150, -140, 60], {'A', 'B', 'D', 'NP'}, 'FontSize', fs, 'FontWeight', 'b','Color', c, 'HorizontalAlignment', 'center')
set(gca, 'FontSize', fs)
set(ax2, 'Color', 'none', 'Box', 'off', 'XColor', 'none', 'YColor', 'none');


%% Band passing to isolate 3 - 20 day variability

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
slp_grid.highpassed = NaN(size(slp_grid.pressure_anom));
slp_grid.lowpassed = NaN(size(slp_grid.pressure_anom));
slp_grid.bandpassed = NaN(size(slp_grid.pressure_anom));
for i = 1:size(slp_grid.pressure_anom, 1)
    for j = 1:size(slp_grid.pressure_anom, 2)
        [slp_grid.highpassed(i,j,:), ~] = highpassing_butterworth(slp_grid.pressure_anom(i,j,:), sampling_interval, cutoff_lower);
        [tmp, slp_grid.lowpassed(i,j,:)] = highpassing_butterworth(slp_grid.pressure_anom(i,j,:), sampling_interval, cutoff_upper);
        [~, slp_grid.bandpassed(i,j,:)] = highpassing_butterworth(tmp, sampling_interval, cutoff_lower);
    end
end

%% Plotting correlation coefficient between BPR and Arctic SLP (maps)

mooring_no = 'np';

fig = figure('Position', [10 10 1200 500]);
tiledlayout(1,3,'TileSpacing', 'none', 'Padding', 'none')
fs = 12; lw = 2;
sgtitle('Relationship between OBP at Mooring ' + string(upper(mooring_no)) + ' and Arctic SLP')

gridded = slp_grid;
for m = 1:3
    if m == 1
        X = gridded.highpassed; % gridded time series
        Y = bpr.(string(mooring_no) + '_highpassed'); % geographically isolated time series
        t = slp.datetime; % shared time grid
        title_str = '< 3 days';
    elseif m == 2
        X = gridded.bandpassed; % gridded time series
        Y = bpr.(string(mooring_no) + '_bandpassed'); % geographically isolated time series
        t = slp.datetime; % shared time grid
        title_str = '3 - 20 days';
    elseif m == 3
        X = gridded.lowpassed; % gridded time series
        Y = bpr.(string(mooring_no) + '_lowpassed'); % geographically isolated time series
        t = slp.datetime; % shared time grid
        title_str = '> 20 days';
    end

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

            % ind = starts(1):ends(end);
            % x = x(ind);
            % y = y(ind);

            %%% Computing significance levels
            if sum(isnan(x)) > 0.1*length(x)
                continue
            end
            % x(isnan(x)) = 0;
            % y(isnan(y)) = 0;
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
    R(flagit) = NaN;
    

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
    clim(0.5.*[-1 1])
    colormap(cmocean('balance'))
    [~, c] = colornames('Crayola', 'white');
    textm([75, 78, 74, 90], [-150, -150, -140, 60], {'A', 'B', 'D', 'NP'}, 'FontSize', 8, 'FontWeight', 'b','Color', c)
    title(string(title_str))

    set(gca, 'FontSize', fs)
    set(ax, 'Color', 'none', 'Box', 'off', 'XColor', 'none', 'YColor', 'none');
end
