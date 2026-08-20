%% Formatting and detiding BPR data each mooring site
clear; clc;

%%% Mooring location
mooring_no = 'd';

if strcmp(mooring_no, 'a')
    lat = 75;
    lon = -150;
elseif strcmp(mooring_no, 'b')
    lat = 78;
    lon = -150;
elseif strcmp(mooring_no, 'd')
    lat = 74;
    lon = -140;
end

%%% Creating time grid
dt = minutes(30);
first_year = 2003;
last_year = 2025;
tgrid = datetime(first_year, 1, 1, "TimeZone", "UTC"):dt:datetime(last_year+1, 1, 1, "TimeZone", "UTC");
tgrid = tgrid(1:end-1);
tgrid = tgrid(:);

%%% Year-pairings for BPR files
yrs = {'03','04';'04','05';'05','06'; '06','07';'07','08'; '08','09'; '09','10'; '10','11';'11','12';'12','13';'13','14';'14','15';'15','16';'16','17'; '17','18'; '18','21'; '21', '22'; '22', '23'; '23', '24'; '24', '25'};

%%% Getting bpr files
u = 1;
for i = 1:size(yrs, 1)
   
    bpr_files = dir('C:\Users\jak279\OneDrive - Yale University\Research\Data\mooring_data\mooring_' + string(mooring_no) + '\bg' + string(yrs{i,1}) + string(yrs{i,2}) + '*.dat');
    if ~isempty(bpr_files)

        bpr_tmp = readtable(string(bpr_files.folder) + '\' + string(bpr_files.name));
        bpr_tmp = renamevars(bpr_tmp, ["pressure_dbar_", "temperature_C_"], ["pressure", "temperature"]);
        bpr_tmp.datetime = datetime(string(num2str(bpr_tmp.x_date)) + ' ' + string(num2str(bpr_tmp.time_UTC_, '%04.f')), 'InputFormat',"yyyyMMdd HHmm", "TimeZone", "UTC");

        %%% Formatting to column vector
        bpr_final.datetime = bpr_tmp.datetime(:);
        bpr_final.pressure = bpr_tmp.pressure(:);

        %%% Removing bad data from the 2005-2006 deployment at Mooring D
        if strcmp(mooring_no, 'd') & strcmp(yrs{i,1}, '05')
            ind = bpr_final.datetime < datetime(2005,9,10,'TimeZone','UTC') | bpr_final.datetime > datetime(2006,5,15,14,0,0,'TimeZone', 'UTC');
            bpr_final.pressure(ind) = NaN;
        end

        %%% Computing BPR anomalies for each deployment
        bpr_final.pressure_anomaly = bpr_final.pressure - mean(bpr_final.pressure, 'omitnan');

        %%% Dividing the 2018-2021 deployment into ~1 year segments
        if strcmp(yrs{i,1}, '18')
            n = floor(length(bpr_final.datetime) / 3);
            dividing_indices = [0, n, 2*n, length(bpr_final.datetime)];

            for m = 1:3
                idx = (dividing_indices(m)+1):dividing_indices(m+1);
                bpr_final_seg.pressure = bpr_final.pressure(idx);
                bpr_final_seg.pressure_anomaly = bpr_final_seg.pressure - mean(bpr_final_seg.pressure, 'omitnan');
                bpr_final_seg.datetime = bpr_final.datetime(idx);

                bpr_all(u,1) = bpr_final_seg; 
                u = u + 1;
            end
            clear bpr_final_seg idx n m dividing_indices

            continue
        end
    
        %%% Saving data to a structure
        bpr_all(u,1) = bpr_final;
        clear bpr_tmp bpr_final

    end

    u = u + 1;
    clear bpr_files
end

%%% Extracting first and last day of deployment
for i = 1:size(bpr_all, 1)
    t = bpr_all(i).datetime(isfinite(bpr_all(i).datetime));
    bpr_all(i,1).start = min(t);
    bpr_all(i,1).end = max(t);
end


%%% Interpolating to get BPR data on time grid
clear bpr_interp_all t_arr
for ii = 1:size(bpr_all,1)

    %%% Skipping if empty
    if isempty(bpr_all(ii,1).pressure)
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
    bpr_interp.datetime = t_arr.datetime;
    bpr_interp.pressure = NaN(1, length(t_arr.datetime));
    bpr_interp.pressure_anomaly = NaN(1, length(t_arr.datetime));

    %%% Looping over continuous sections of data
    good = isfinite(bpr.pressure);
    starts = find(diff([0; good]) == 1);
    ends   = find(diff([good; 0]) == -1);
    for uu = 1:length(starts)
        ind = starts(uu):ends(uu);
        idx = bpr_interp.datetime >= bpr.datetime(starts(uu)) & bpr_interp.datetime <= bpr.datetime(ends(uu));
        bpr_interp.pressure(idx) = interp1(bpr.datetime(ind),  bpr.pressure(ind), bpr_interp.datetime(idx));
        bpr_interp.pressure_anomaly(idx) = interp1(bpr.datetime(ind),  bpr.pressure_anomaly(ind), bpr_interp.datetime(idx));
    end
    bpr_interp.start = bpr_all(ii,1).start;
    bpr_interp.end = bpr_all(ii,1).end;

    %%% Removing NaNs
    idx = find(~isnan(bpr_interp.pressure_anomaly), 1, 'first');
    idy = find(~isnan(bpr_interp.pressure_anomaly), 1, 'last');
    bpr_interp.datetime = bpr_interp.datetime(idx:idy);
    bpr_interp.pressure = bpr_interp.pressure(idx:idy);
    bpr_interp.pressure_anomaly = bpr_interp.pressure_anomaly(idx:idy);

    %%% Formatting into column vector
    bpr_interp.pressure = bpr_interp.pressure(:);
    bpr_interp.pressure_anomaly = bpr_interp.pressure_anomaly(:);

    %%% Saving
    bpr_interp_all(ii,1) = bpr_interp;

end

clear bpr bpr_interp ii i t1 t2 t_arr ind

%%% Renaming data
bpr_highres_all = bpr_interp_all;
clear bpr_interp_all


%%

%%% Detiding annual segments using t-tide
for i = 1:length(bpr_highres_all)
    if ~isempty(bpr_highres_all(i).pressure_anomaly)


        t = bpr_highres_all(i).datetime;
        y = bpr_highres_all(i).pressure_anomaly;
        ttide_dt = (datenum(t(2)) - datenum(t(1))) * 24; % units of hours

        %%% Finding continuous segments
        good = isfinite(y);
        starts = find(diff([0; good]) == 1);
        ends   = find(diff([good; 0]) == -1);

        %%% Detiding via t-tide
        detided = NaN(size(y));
        for u = 1:length(starts)
            ind = starts(u):ends(u);

            [nameu,fu,tidecon,pout]=t_tide(y(ind),...
                'interval', ttide_dt, ... % sampling interval
                'start', datenum(t(starts(u))),...
                'latitude', lat); % latitude of observation);

            detided(ind) = y(ind) - pout;

        end

        bpr_highres_all(i).detided = detided;

    end
end


%%% Combining high-res data into one long time series
clear combined
bpr = bpr_highres_all;
combined.datetime = bpr(1).datetime;
combined.pressure = bpr(1).pressure;
combined.pressure_anomaly = bpr(1).pressure_anomaly;
combined.detided = bpr(1).detided;
if size(bpr_highres_all,1) > 1
    for i = 2:size(bpr_highres_all, 1)
        combined.datetime = [combined.datetime; bpr(i).datetime];
        combined.pressure = [combined.pressure; bpr(i).pressure];
        combined.pressure_anomaly = [combined.pressure_anomaly; bpr(i).pressure_anomaly];
        combined.detided = [combined.detided; bpr(i).detided];
    end
end
clear i

%%% Add NaNs to fill data gaps
combined_final.datetime = tgrid;
combined_final.pressure = NaN(size(tgrid));
combined_final.pressure_anomaly = NaN(size(tgrid));
combined_final.detided = NaN(size(tgrid));
for i = 1:length(combined_final.datetime)
    [~,idx] = ismember(combined_final.datetime(i), combined.datetime);
    if idx == 0
        continue
    end
    combined_final.pressure(i) = combined.pressure(idx);
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
clear dt_string bpr bpr_all 

%% Putting BPR data (raw + detided) on coarser time grid

%%% Creating coarser time grid
dt = hours(6);
tgrid = datetime(first_year, 1, 1, "TimeZone", "UTC"):dt:datetime(last_year+1, 1, 1, "TimeZone", "UTC");
if dt == days(7)
    tgrid = datetime(2003, 1, 5, "TimeZone", "UTC"):dt:datetime(2015, 1, 12, "TimeZone", "UTC");
end
tgrid = tgrid(1:end-1);
tgrid = tgrid(:);

%%% Running mean 
clear bpr_all t_arr
for ii = 1:size(bpr_highres_all,1)

    %%% Skipping if empty
    if isempty(bpr_highres_all(ii,1).pressure)
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
    [bpr_interp.datetime, bpr_interp.pressure] = running_mean(bpr.datetime, bpr.pressure, t_arr.datetime);
    [bpr_interp.datetime, bpr_interp.pressure_anomaly] = running_mean(bpr.datetime, bpr.pressure_anomaly, t_arr.datetime);
    [bpr_interp.datetime, bpr_interp.detided] = running_mean(bpr.datetime, bpr.detided, t_arr.datetime);

    % bpr_interp.pressure = NaN(length(t_arr.datetime),1);
    % bpr_interp.pressure_anomaly = NaN(length(t_arr.datetime),1);
    % bpr_interp.detided = NaN(length(t_arr.datetime),1);
    % for i = 2:length(t_arr.datetime) - 1
    %     ind = (datenum(bpr.datetime) > (datenum(t_arr.datetime(i)) - (dt/2))) & (datenum(bpr.datetime) <= (datenum(t_arr.datetime(i)) + (dt/2)));
    %     bpr_interp.pressure(i) = mean(bpr.pressure(ind));
    %     bpr_interp.pressure_anomaly(i) = mean(bpr.pressure_anomaly(ind));
    %     bpr_interp.detided(i) = mean(bpr.detided(ind));
    % end
    % bpr_interp.datetime = t_arr.datetime;
    bpr_interp.start = bpr.start;
    bpr_interp.end = bpr.end;

    %%% Removing NaNs
    idx = ~isnan(bpr_interp.detided);
    bpr_interp.datetime = bpr_interp.datetime(idx);
    bpr_interp.pressure = bpr_interp.pressure(idx);
    bpr_interp.pressure_anomaly = bpr_interp.pressure_anomaly(idx);
    bpr_interp.detided = bpr_interp.detided(idx);

    %%% Saving data
    bpr_all(ii,1) = bpr_interp;

end

%%% Combining daily data into one long time series
clear combined
bpr = bpr_all;
combined.datetime = bpr(1).datetime;
combined.pressure = bpr(1).pressure;
combined.pressure_anomaly = bpr(1).pressure_anomaly;
combined.detided = bpr(1).detided;
if size(bpr_all,1) > 1
    for i = 2:size(bpr_all, 1)
        combined.datetime = [combined.datetime; bpr(i).datetime];
        combined.pressure = [combined.pressure; bpr(i).pressure];
        combined.pressure_anomaly = [combined.pressure_anomaly; bpr(i).pressure_anomaly];
        combined.detided = [combined.detided; bpr(i).detided];
    end
end
clear i

%%% Add NaNs to fill data gaps
combined_final.datetime = tgrid;
combined_final.pressure = NaN(size(tgrid));
combined_final.pressure_anomaly = NaN(size(tgrid));
combined_final.detided = NaN(size(tgrid));
for i = 1:length(combined_final.datetime)
    [~,idx] = ismember(combined_final.datetime(i), combined.datetime);
    if idx == 0
        continue
    end
    combined_final.pressure(i) = combined.pressure(idx);
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
end

clear ii i t1 t2 t_arr ind

%%% Saving data
bpr = combined_final;
bpr.lat = lat;
bpr.lon = lon;
dt_string = char(dt);
dt_string = string(dt_string(~isspace(dt_string)));
save('C:\Users\jak279\OneDrive - Yale University\Research\Data\OBP_interp\BPR_mooring_' + string(mooring_no) + '_' + string(dt_string) + '.mat', 'bpr', 'bpr_all')

%%

figure('Position', [10 10 1000 700])
tiledlayout(2,2)
i = 3;
lw = 2;
sgtitle('Mooring ' + string(upper(mooring_no)))

nexttile([1 2])
plot(bpr_highres_all(i).datetime, bpr_highres_all(i).pressure_anomaly, 'Color', 'k', 'LineWidth', lw)
hold on
plot(bpr_highres_all(i).datetime, bpr_highres_all(i).detided, 'Color', [0.5 0.5 0.5], 'LineWidth', lw)
grid on
xlim([bpr_highres_all(i).start bpr_highres_all(i).end])

nexttile()
spectral_estimate_pwelch(bpr_highres_all(i).datetime, bpr_highres_all(i).pressure_anomaly, 'Color', 'k', 'MakeFigure', 1)
spectral_estimate_pwelch(bpr_highres_all(i).datetime, bpr_highres_all(i).detided, 'Color', [0.5 0.5 0.5], 'MakeFigure', 1)
title('Pwelch')

P = 10;
nexttile()
spectral_estimate_jlab(bpr_highres_all(i).datetime, bpr_highres_all(i).pressure_anomaly, 'P', P, 'Color', 'k', 'MakeFigure', 1)
spectral_estimate_jlab(bpr_highres_all(i).datetime, bpr_highres_all(i).detided, 'P', P, 'Color', [0.5 0.5 0.5], 'MakeFigure', 1)
title('Slepian Multi-Taper')