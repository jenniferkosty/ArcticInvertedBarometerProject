%% Script to format BPR time series from PSMSL

clc; clear;

%%

mooring_no = 'northpacific';

if strcmp(mooring_no, 'northatlantic')
    directory = 'C:\Users\jak279\OneDrive - Yale University\Research\Data\PSMSL_BPR\DART_44401_hrp\data\';
    lat = 37.551;
    lon = -49.985;
elseif strcmp(mooring_no, 'northpacific')
    directory = 'C:\Users\jak279\OneDrive - Yale University\Research\Data\PSMSL_BPR\NDBC_46402_hrp\data\';
    lat = 50.8881;
    lon = -164.3138;
elseif strcmp(mooring_no, 'hawaii')
    directory = 'C:\Users\jak279\OneDrive - Yale University\Research\Data\PSMSL_BPR\NDBC_51407_hrp\data\';
    lat = 19.5777;
    lon = -156.5393;
elseif strcmp(mooring_no, 'drake')
    directory = 'C:\Users\jak279\OneDrive - Yale University\Research\Data\PSMSL_BPR\DPN_DEEP_hrp\data\';
    lat = -56.035;
    lon = -57.9658;
end

%%% dt for coarse time grid
dt_coarse = days(1);

file_info = dir(fullfile(directory, '*.txt'));
u = 1;
for m = 1:length(file_info)
    fid = fopen(directory + string(file_info(m).name), 'r');

    % Skip header lines (adjust number as needed)
    for k = 1:20
        fgetl(fid);
    end

    % Read the 7 columns: int, int, int, int, float, float, int
    data = textscan(fid, '%d %d %d %d %f %f %f %f %f');
    fclose(fid);

    % Extract into named variables
    Recno    = data{1};   % Record number
    Fl       = data{2};   % Flag (0=good, 1=bad/missing)
    Year     = data{3};   % Year
    Day      = data{4};   % Day of year
    Hour     = data{5};   % Hour of day
    Pres     = data{6};   % Pressure [mbar]
    ResDrDft = data{9};   % Residual pressure (tidal + drift removal) [mbar]

    % Converting to datetime
    bpr_tmp.datetime = datetime(Year, 1, 1, 'TimeZone', 'UTC') + days(Day - 1) + hours(Hour);
    bpr_tmp.detided = ResDrDft ./ 100; %%% mb to dbar
    bpr_tmp.detided(Fl == 1) = NaN;

    %%% Getting approximately 1 year segments
    total_years = years(bpr_tmp.datetime(end) - bpr_tmp.datetime(1));
    if total_years > 1.1
        num_segments = round(total_years);

        % Generate evenly spaced breakpoints
        breakpoints = bpr_tmp.datetime(1) + (0:num_segments) * (bpr_tmp.datetime(end) - bpr_tmp.datetime(1)) / num_segments;

        % Build start and end date arrays
        start_dates = breakpoints(1:end-1);
        end_dates   = breakpoints(2:end);
    else
        start_dates = bpr_tmp.datetime(1);
        end_dates = bpr_tmp.datetime(end);
    end

    %%% Creating ~1-year segments of bpr data
    clear bpr_all
    for i = 1:length(start_dates)
        ind = bpr_tmp.datetime >= start_dates(i) & bpr_tmp.datetime <= end_dates(i);
        bpr_highres_all(u,1).datetime = bpr_tmp.datetime(ind);
        bpr_highres_all(u,1).detided = bpr_tmp.detided(ind) - mean(bpr_tmp.detided(ind), 'omitnan');
        bpr_highres_all(u,1).start = start_dates(i);
        bpr_highres_all(u,1).end = end_dates(i);
        u = u + 1;
    end
end

%%% Creating time grid
dt = hours(1);
first_year = 2003;
last_year = 2025;
tgrid = datetime(first_year, 1, 1, "TimeZone", "UTC"):dt:datetime(last_year+1, 1, 1, "TimeZone", "UTC");
tgrid = tgrid(1:end-1);
tgrid = tgrid(:);

%%% Combining segments
bpr = bpr_highres_all;

combined.datetime = bpr(1).datetime;
combined.detided  = bpr(1).detided;
if size(bpr,1) > 1
    for i = 2:size(bpr,1)
        combined.datetime = [combined.datetime; bpr(i).datetime];
        combined.detided  = [combined.detided;  bpr(i).detided];
    end
end

%%% Add NaNs to fill data gaps
combined_final.datetime = tgrid;
combined_final.detided = NaN(size(tgrid));
for i = 1:length(combined_final.datetime)
    [~,idx] = ismember(combined_final.datetime(i), combined.datetime);
    if idx == 0
        continue
    end
    combined_final.detided(i) = combined.detided(idx);
end

%%% Adding column with start/end dates
combined_final.start = NaT(1, length(bpr), 'TimeZone', 'UTC');
combined_final.end = NaT(1, length(bpr), 'TimeZone', 'UTC');
for i = 1:length(bpr)
    if isempty(bpr(i).start)
        continue
    end
    combined_final.start(i) = bpr(i).start;
    combined_final.end(i) = bpr(i).end;
end

bpr_highres = combined_final;
bpr_highres.lat = lat;
bpr_highres.lon = lon;
clear combined_final idx

%%% Saving data
dt_string = char(dt);
dt_string = string(dt_string(~isspace(dt_string)));
bpr = bpr_highres;
bpr_all = bpr_highres_all;
save('C:\Users\jak279\OneDrive - Yale University\Research\Data\OBP_interp\BPR_mooring_' + string(mooring_no) + '_' + string(dt_string) + '.mat', 'bpr', 'bpr_all')

%%% Putting BPR data on coarser time grid

%%% Creating coarser time grid
dt = dt_coarse;
tgrid = datetime(first_year, 1, 1, "TimeZone", "UTC"):dt:datetime(last_year+1, 1, 1, "TimeZone", "UTC");
tgrid = tgrid(1:end-1);
tgrid = tgrid(:);

%%% Running mean 
clear bpr_all t_arr
for ii = 1:size(bpr_highres_all,1)

    %%% Skipping if empty
    if isempty(bpr_highres_all(ii,1).detided)
        continue
    end

    %%% Getting time array for each year
    t1 = datetime(year(bpr_highres_all(ii).start), 1, 1, "TimeZone", "UTC");
    t2 = datetime(year(bpr_highres_all(ii).end)+1, 12, 31, "TimeZone", "UTC");
    t_arr.datetime = tgrid(tgrid >= t1 & tgrid <= t2);

    %%% Selecting bpr data for 1 year
    bpr = bpr_highres_all(ii,1);

    %%% Getting bpr data on time grid
    clear bpr_interp
    bpr_interp.datetime = t_arr.datetime;
    bpr_interp.detided = NaN(length(t_arr.datetime),1);
    for i = 2:length(t_arr.datetime) - 1
        ind = (datenum(bpr.datetime) > (datenum(t_arr.datetime(i)) - (dt/2))) & (datenum(bpr.datetime) <= (datenum(t_arr.datetime(i)) + (dt/2)));
        bpr_interp.detided(i) = mean(bpr.detided(ind));
    end
    bpr_interp.start = bpr.start;
    bpr_interp.end = bpr.end;

    %%% Removing NaNs
    idx = ~isnan(bpr_interp.detided);
    bpr_interp.datetime = bpr_interp.datetime(idx);
    bpr_interp.detided = bpr_interp.detided(idx);
    bpr_all(ii,1) = bpr_interp;

end

%%% Combining coarse data into one long time series
clear combined
bpr = bpr_all;
combined.datetime = bpr(1).datetime;
combined.detided = bpr(1).detided;
if size(bpr_all,1) > 1
    for i = 2:size(bpr_all, 1)
        combined.detided = [combined.detided; bpr(i).detided];
        combined.datetime = [combined.datetime; bpr(i).datetime];
    end
end
clear i

%%% Add NaNs to fill data gaps
combined_final.datetime = tgrid;
combined_final.detided = NaN(size(tgrid));
for i = 1:length(combined_final.datetime)
    [~,idx] = ismember(combined_final.datetime(i), combined.datetime);
    if idx == 0
        continue
    end
    combined_final.detided(i) = combined.detided(idx);
end

%%% Adding column with start/end dates
combined_final.start = NaT(1, length(bpr), 'TimeZone', 'UTC');
combined_final.end = NaT(1, length(bpr), 'TimeZone', 'UTC');
for i = 1:length(bpr)
    if isempty(bpr(i).start)
        continue
    end
    combined_final.start(i) = bpr(i).start;
    combined_final.end(i) = bpr(i).end;
end

clear ii i t1 t2 t_arr ind

%%% Saving data
bpr = combined_final;
bpr.lat = lat;
bpr.lon = lon;
dt_string = char(dt);
dt_string = string(dt_string(~isspace(dt_string)));
save('C:\Users\jak279\OneDrive - Yale University\Research\Data\OBP_interp\BPR_mooring_' + string(mooring_no) + '_' + string(dt_string) + '.mat', 'bpr', 'bpr_all');

%% Map figure showing mooring location

figure()
for i = 1:3
    if i == 1
        lat(i) = 37.551;
        lon(i) = -49.985;
    elseif i == 2
        lat(i) = 50.8881;
        lon(i) = -164.3138;
    elseif i == 3
        lat(i) = 19.5777;
        lon(i) = -156.5393;
    end
    geoscatter(lat, lon, 'k', 'filled', '')
end
geobasemap colorterrain