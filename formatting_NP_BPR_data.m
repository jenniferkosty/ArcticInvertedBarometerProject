%%% Script to format NP BPR record
clc; clear;

%%

%%% Selecting version
version = 2;

%%% Setting coarse grid dt
dt_coarse = hours(6);

%%% Latitude of record
lat = 89.25;
lon = 60.35;

%%% First section: 2005-2010
path = "C:\Users\jak279\OneDrive - Yale University\Research\Data\NP_mooring\data\North_Pole_5-Year_Arctic_Bottom_Pressure_2005-2010\ABPR_NorthPoleDATA_5years.txt";
data = readtable(path);
time_a = datetime(data.(1), data.(2), data.(3), data.(4), data.(5), 0, 'TimeZone', 'UTC');
datenum_a = datenum(time_a);
bpr_a = data.(9); % bottom pressure anomalies relative to the total mean = raw - drift - averaged pressure.
%bpr_a = bpr_a ./ 100; % meters of equivalent sea water

%%% These lines convert to dbar!
bpr_a = bpr_a + 430094.629; % average depth in cm;
bpr_a = bpr_a ./ 100; % converting to meters
bpr_a = gsw_p_from_z(-bpr_a, lat);
bpr_a = bpr_a - mean(bpr_a, 'omitnan');

% %%% First section: 2005-2010 (hourly time grid) - USED FOR COMPARISON
% path = "C:\Users\jak279\OneDrive - Yale University\Research\Data\NP_mooring\data\North_Pole_5-Year_Arctic_Bottom_Pressure_2005-2010\ABPR_NorthPoleDATA_5years_HourlyData.txt";
% data = readtable(path);
% time_a_hourly = datetime(data.(1), data.(2), data.(3), data.(4), data.(5), 0, 'TimeZone', 'UTC');
% datenum_a_hourly = datenum(time_a_hourly);
% tides_a_hourly = data.(8); % tidal signal using matlab program t_tide.m 
% bpr_a_hourly_detided = data.(9); % bottom pressure anomalies de-tided and de-meaned.

%%% Second section: 2010-2014
path = "C:\Users\jak279\OneDrive - Yale University\Research\Data\NP_mooring\data\North_Pole_Arctic_Bottom_Pressure_2010-2014\ABPR5_2010to2014data.txt";
data = readtable(path);
time_b = datetime(data.(1), data.(2), data.(3), data.(4), data.(5), data.(6), 'TimeZone', 'UTC');
datenum_b = datenum(time_b);
bpr_b = data.(9); % de-drifted ocean bottom pressure anomalies relative to the total mean = raw - drift - averaged pressure
tides_b = data.(10); % tidal amplitude (obtained from the matlab program T_tide, by Pawlowicz, 2002)
bpr_b_detided = data.(11); % de-tided, de-drifted, ocean bottom pressure anomalies.
%bpr_b = bpr_b ./ 100; % meters of equivalent sea water

%%% These lines convert to dbar!
bpr_b = bpr_b + 419710.0; % average depth in cm;
bpr_b = bpr_b ./ 100; % converting to meters
bpr_b = gsw_p_from_z(-bpr_b, lat);
bpr_b = bpr_b - mean(bpr_b, 'omitnan');

clear path data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Version 1 of Processing NP data %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if version == 1

% Just process the 2 multi-year deployments

%%% Hourly time grid
dt = hours(1);
first_year = 2003;
last_year = 2025;
tgrid = datetime(first_year, 1, 1, "TimeZone", "UTC"):dt:datetime(last_year+1, 1, 1, "TimeZone", "UTC");
tgrid = tgrid(1:end-1);
tgrid = tgrid(:);

%%% Interpolating onto hourly grid
bpr_all = [bpr_a; bpr_b];
time_all = [time_a; time_b];
bpr_highres.datetime = tgrid;
bpr_highres.pressure_anomaly = interp1(time_all, bpr_all, tgrid);
ind = bpr_highres.datetime > time_a(end) & bpr_highres.datetime < time_b(1);
bpr_highres.pressure_anomaly(ind) = NaN;
bpr_highres.lat = lat;
bpr_highres.lon = lon;

%%% Detiding
unique_years = unique(year(bpr_highres.datetime));
bpr_highres.detided = NaN(size(bpr_highres.pressure_anomaly));
for i = 1:length(unique_years)
    ind = find(year(bpr_highres.datetime) == unique_years(i));
    x = bpr_highres.datetime(ind);
    y = bpr_highres.pressure_anomaly(ind);
    dt_ttide = (datenum(x(2)) - datenum(x(1))) * 24; % converting to units of hours
    if sum(isnan(y)) > 0.5*length(y)
        continue
    end
    [~,pout]=t_tide(y,...
        'interval', dt_ttide, ... % sampling interval
        'start', datenum(x(1)),...
        'latitude', lat); % latitude of observation);
    bpr_highres.detided(ind) = bpr_highres.pressure_anomaly(ind) - pout;
end
bpr_highres.start = [time_a(1); time_b(1)];
bpr_highres.end = [time_a(end); time_b(end)];
bpr = bpr_highres;

%%% Saving highres data
dt_string = char(dt);
dt_string = string(dt_string(~isspace(dt_string)));
save('C:\Users\jak279\OneDrive - Yale University\Research\Data\OBP_interp\BPR_mooring_np_' + string(dt_string) + '.mat', 'bpr')

%%% Putting North Pole data on a coarser grid

%%% Creating 6-hourly time grid
dt = dt_coarse;
tgrid = datetime(first_year, 1, 1, "TimeZone", "UTC"):dt:datetime(last_year+1, 1, 1, "TimeZone", "UTC");
tgrid = tgrid(1:end-1);
tgrid = tgrid(:);

%%% Running mean to get bpr data on coarser grid
clear bpr_interp_hourly
bpr_interp.pressure_anomaly = NaN(size(tgrid));
bpr_interp.detided = NaN(size(tgrid));
for i = 2:length(tgrid) - 1
    ind = (bpr_highres.datetime >= (tgrid(i) - (dt/2))) & (bpr_highres.datetime < (tgrid(i) + (dt/2)));
    bpr_interp.pressure_anomaly(i) = mean(bpr_highres.pressure_anomaly(ind), 'omitnan');
    bpr_interp.detided(i) = mean(bpr_highres.detided(ind), 'omitnan');
end
bpr_interp.datetime = tgrid;
bpr_interp.start = bpr_highres.start;
bpr_interp.end = bpr_highres.end;
bpr = bpr_interp;

dt_string = char(dt);
dt_string = string(dt_string(~isspace(dt_string)));
save('C:\Users\jak279\OneDrive - Yale University\Research\Data\OBP_interp\BPR_mooring_np_' + string(dt_string) + '.mat', 'bpr')

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Version 2 of Processing NP data %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if version == 2

% Treat data exactly like BG BPR data (compute anomalies on ~year long
% deployments)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% OPTION A: Use exact start/end points from BG BPR deployments %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%% Loading BG BPR to get start/end points
% moorings_all = {'a', 'b', 'd'};
% for i = 1:length(moorings_all)
% 
%     %%% Selecting mooring
%     mooring_no = moorings_all{i};
% 
%     %%% Loading BPR data
%     data = load('C:\Users\jak279\OneDrive - Yale University\Research\Data\OBP_interp\BPR_mooring_' + string(mooring_no) + '_6hr.mat'); % BPR
%     bpr.(string(mooring_no) + '_start') = data.bpr.start;
%     bpr.(string(mooring_no) + '_end') = data.bpr.end;
% end
% 
% start_dates = mean([bpr.a_start; bpr.b_start; bpr.d_start], 1, 'omitnan');
% end_dates = mean([bpr.a_end; bpr.b_end; bpr.d_end], 1, 'omitnan');
% clear bpr data mooring_no moorings_all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% OPTION B: Use ~1 year long segments, maximizing data usage %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

time_tmp = [time_a; time_b];
total_years = years(time_tmp(end) - time_tmp(1));
if total_years > 1.5
    num_segments = round(total_years);

    % Generate evenly spaced breakpoints
    breakpoints = time_tmp(1) + (0:num_segments) * (time_tmp(end) - time_tmp(1)) / num_segments;

    % Build start and end date arrays
    start_dates = breakpoints(1:end-1);
    end_dates   = breakpoints(2:end);
end
clear time_tmp total_years num_segments breakpoints

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% BOTH OPTIONS PROCEED THE SAME FROM HERE %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Isolating NP data data onto 1 year segments
clear bpr_all
for i = 1:length(start_dates)
    if end_dates(i) < time_a(end)
        ind = time_a >= start_dates(i) & time_a <= end_dates(i);
        bpr_all(i,1).datetime = time_a(ind);
        bpr_all(i,1).pressure_anomaly = bpr_a(ind) - mean(bpr_a(ind), 'omitnan');
        if sum(ind) ~= 0 
            bpr_all(i,1).duration = (bpr_all(i).datetime(end) - bpr_all(i).datetime(1)) / (end_dates(i) - start_dates(i)) * 100;
        else
            bpr_all(i,1).duration = 0;
        end
    else
        ind = time_b >= start_dates(i) & time_b <= end_dates(i);
        bpr_all(i,1).datetime = time_b(ind);
        bpr_all(i,1).pressure_anomaly = bpr_b(ind) - mean(bpr_b(ind), 'omitnan');
        if sum(ind) ~= 0
            bpr_all(i,1).duration = (bpr_all(i).datetime(end) - bpr_all(i).datetime(1)) / (end_dates(i) - start_dates(i)) * 100;
        else
            bpr_all(i,1).duration = 0;
        end
    end
    bpr_all(i,1).start = start_dates(i);
    bpr_all(i,1).end = end_dates(i);
end

%%% Creating time grid
dt = hours(1);
first_year = 2003;
last_year = 2025;
tgrid = datetime(first_year, 1, 1, "TimeZone", "UTC"):dt:datetime(last_year+1, 1, 1, "TimeZone", "UTC");
tgrid = tgrid(1:end-1);
tgrid = tgrid(:);

%%% Interpolating to get BPR data on time grid
clear bpr_interp_all
for ii = 1:size(bpr_all,1)

    %%% Skipping if empty
    if isempty(bpr_all(ii).pressure_anomaly)
        bpr_interp_all(ii,1) = bpr_all(ii);
        continue
    end

    %%% Getting time segment for each year
    t1 = datetime(year(bpr_all(ii).start), 1, 1, "TimeZone", "UTC");
    t2 = datetime(year(bpr_all(ii).end)+1, 1, 1, "TimeZone", "UTC");
    t_arr.datetime = tgrid(tgrid >= t1 & tgrid < t2);

    %%% Selecting bpr data for 1 year
    bpr = bpr_all(ii,1);

    %%% Interpolating BPR data
    clear bpr_interp
    bpr_interp.pressure_anomaly = NaN(size(t_arr.datetime));
    bpr_interp.datetime = t_arr.datetime;

    %%% Looping over continuous sections of data
    good = isfinite(bpr.pressure_anomaly);
    starts = find(diff([0; good]) == 1);
    ends   = find(diff([good; 0]) == -1);
    for uu = 1:length(starts)
        ind = starts(uu):ends(uu);
        idx = bpr_interp.datetime >= bpr.datetime(starts(uu)) & bpr_interp.datetime <= bpr.datetime(ends(uu));
        bpr_interp.pressure_anomaly(idx) = interp1(bpr.datetime(ind),  bpr.pressure_anomaly(ind), bpr_interp.datetime(idx));
    end
    bpr_interp.start = bpr_all(ii,1).start;
    bpr_interp.end = bpr_all(ii,1).end;
    bpr_interp.duration = bpr_all(ii,1).duration;

    %%% Removing NaNs
    idx = find(~isnan(bpr_interp.pressure_anomaly), 1, 'first');
    idy = find(~isnan(bpr_interp.pressure_anomaly), 1, 'last');
    bpr_interp.pressure_anomaly = bpr_interp.pressure_anomaly(idx:idy);
    bpr_interp.datetime = bpr_interp.datetime(idx:idy);

    %%% Formatting into column vector
    bpr_interp.pressure_anomaly = bpr_interp.pressure_anomaly(:);

    %%% Saving
    bpr_interp_all(ii,1) = bpr_interp;

end

clear bpr bpr_interp ii i t1 t2 t_arr ind idx idy good starts ends uu
%%
%%% Renaming data
bpr_highres_all = bpr_interp_all;
clear bpr_interp_all
colors = orderedcolors("gem");

%%% Detiding annual segments using t-tide
for i = 1:length(bpr_highres_all)
    if ~isempty(bpr_highres_all(i).pressure_anomaly)

        t = bpr_highres_all(i).datetime;
        y = bpr_highres_all(i).pressure_anomaly;
        t_num = datenum(t);
        ttide_dt = (t_num(2) - t_num(1)) * 24; % units of hours

        %%% Finding NaNs in segment
        nan_mask = ~isfinite(y);
        gap_starts = find(diff([0; nan_mask]) == 1);
        gap_ends   = find(diff([nan_mask; 0]) == -1);
        gap_lengths_hours = (gap_ends - gap_starts + 1) * ttide_dt;

        %%% Skipping too many NaNs
        if sum(gap_lengths_hours) > 0.5 * ((t_num(end) - t_num(1)) * 24)
            detided = NaN(size(y));
            continue
        end


        %%% Detiding
        [~, ~, ~, pout] = t_tide(y, ...
            'interval', ttide_dt, ...
            'start',    t_num(1), ...
            'latitude', lat);
        detided = y - pout;

        %%% Ensure long gap positions stay NaN 
        detided(nan_mask) = NaN;

        %%% Saving
        bpr_highres_all(i).detided = detided;

    end

    clear detided good seg_starts seg_ends short_gaps idx y_filled g nan_mask gap_starts gap_ends gap_length_hours t_num

end

clear nameu fu tidecon pout ttide_dt t y good starts ends detided u ind i

%%% Combining high-res data into one long time series
clear combined
bpr = bpr_highres_all;
combined.datetime = bpr(1).datetime;
combined.pressure_anomaly = bpr(1).pressure_anomaly;
combined.detided = bpr(1).detided;
if size(bpr_highres_all,1) > 1
    for i = 2:size(bpr_highres_all, 1)
        combined.pressure_anomaly = [combined.pressure_anomaly; bpr(i).pressure_anomaly];
        combined.datetime = [combined.datetime; bpr(i).datetime];
        combined.detided = [combined.detided; bpr(i).detided];
    end
end
clear i

%%% Add NaNs to fill data gaps
combined_final.datetime = tgrid;
combined_final.pressure_anomaly = NaN(size(tgrid));
combined_final.detided = NaN(size(tgrid));
for i = 1:length(combined_final.datetime)
    [~,idx] = ismember(combined_final.datetime(i), combined.datetime);
    if idx == 0
        continue
    end
    combined_final.pressure_anomaly(i) = combined.pressure_anomaly(idx);
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
    combined_final.duration(i) = bpr(i).duration;
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
save('C:\Users\jak279\OneDrive - Yale University\Research\Data\OBP_interp\BPR_mooring_np_' + string(dt_string) + '.mat', 'bpr', 'bpr_all')
clear dt_string bpr bpr_all 

%% Putting BPR data (raw + detided) on coarser time grid

%%% Creating coarser time grid
dt = dt_coarse;
tgrid = datetime(first_year, 1, 1, "TimeZone", "UTC"):dt:datetime(last_year+1, 1, 1, "TimeZone", "UTC");
tgrid = tgrid(1:end-1);
tgrid = tgrid(:);

%%% Running mean 
clear bpr_all t_arr
for ii = 1:size(bpr_highres_all,1)

    %%% Skipping if empty
    if isempty(bpr_highres_all(ii,1).pressure_anomaly)
        bpr_all(ii,1) = bpr_highres_all(ii);
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
    % bpr_interp.pressure_anomaly = NaN(length(t_arr.datetime),1);
    % bpr_interp.detided = NaN(length(t_arr.datetime),1);
    % for i = 2:length(t_arr.datetime) - 1
    %     ind = (datenum(bpr.datetime) > (datenum(t_arr.datetime(i)) - (dt/2))) & (datenum(bpr.datetime) <= (datenum(t_arr.datetime(i)) + (dt/2)));
    %     bpr_interp.pressure_anomaly(i) = mean(bpr.pressure_anomaly(ind));
    %     bpr_interp.detided(i) = mean(bpr.detided(ind));
    % end
    % bpr_interp.datetime = t_arr.datetime;
    [bpr_interp.datetime, bpr_interp.pressure_anomaly] = running_mean(bpr.datetime, bpr.pressure_anomaly, t_arr.datetime);
    [bpr_interp.datetime, bpr_interp.detided] = running_mean(bpr.datetime, bpr.detided, t_arr.datetime);
    bpr_interp.start = bpr.start;
    bpr_interp.end = bpr.end;
    bpr_interp.duration = bpr.duration;

    %%% Removing NaNs
    idx = ~isnan(bpr_interp.pressure_anomaly);
    bpr_interp.pressure_anomaly = bpr_interp.pressure_anomaly(idx);
    bpr_interp.detided = bpr_interp.detided(idx);
    bpr_interp.datetime = bpr_interp.datetime(idx);
    bpr_all(ii,1) = bpr_interp;

end

%%% Combining coarse data into one long time series
clear combined
bpr = bpr_all;
combined.datetime = bpr(1).datetime;
combined.pressure_anomaly = bpr(1).pressure_anomaly;
combined.detided = bpr(1).detided;
if size(bpr_all,1) > 1
    for i = 2:size(bpr_all, 1)
        combined.pressure_anomaly = [combined.pressure_anomaly; bpr(i).pressure_anomaly];
        combined.detided = [combined.detided; bpr(i).detided];
        combined.datetime = [combined.datetime; bpr(i).datetime];
    end
end
clear i

%%% Add NaNs to fill data gaps
combined_final.datetime = tgrid;
combined_final.pressure_anomaly = NaN(size(tgrid));
combined_final.detided = NaN(size(tgrid));
for i = 1:length(combined_final.datetime)
    [~,idx] = ismember(combined_final.datetime(i), combined.datetime);
    if idx == 0
        continue
    end
    combined_final.pressure_anomaly(i) = combined.pressure_anomaly(idx);
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
    combined_final.duration(i) = bpr(i).duration;
end

clear ii i t1 t2 t_arr ind

%%% Saving data
bpr = combined_final;
bpr.lat = lat;
bpr.lon = lon;
dt_string = char(dt);
dt_string = string(dt_string(~isspace(dt_string)));
save('C:\Users\jak279\OneDrive - Yale University\Research\Data\OBP_interp\BPR_mooring_np_' + string(dt_string) + '.mat', 'bpr', 'bpr_all')

end

