%% Formatting ULS daily averaged data

clear; clc;

%%% Setting location
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

%%% Loading SLP data (used for p_ULS computation)
data = load('C:\Users\jak279\OneDrive - Yale University\Research\Data\OBP_interp\SLP_mooring_' + string(mooring_no) + '_1day.mat'); % SLP
slp.datetime = data.slp.datetime;
slp.(string(mooring_no)) = NaN(size(data.slp.pressure));
slp_tmp = data.slp.pressure;
slp_tmp = -gsw_z_from_p(slp_tmp, data.slp.actual_lat);
slp.(string(mooring_no)) = slp_tmp;
clear slp_tmp rho0 g data

%%% Creating time grid
dt = days(1);
first_year = 2003;
last_year = 2025;
tgrid = datetime(first_year, 1, 1, "TimeZone", "UTC"):dt:datetime(last_year+1, 1, 1, "TimeZone", "UTC");
tgrid = tgrid(1:end-1);
tgrid = tgrid(:);

%%% Year-pairings to load in data
yrs = {'03','04';'04','05';'05','06'; '06','07';'07','08'; '08','09'; '09','10'; '10','11';'11','12';'12','13';'13','14';'14','15';'15','16';'16','17'; '17','18'; '18','21'; '21', '22'; '22', '23'; '23', '24'; '24', '25'};

%% Manually removing eddies from ULS data (only need to do once!)

%%% Selecting year
ii = 6;

%%% Loading data
files = dir('C:\Users\jak279\OneDrive - Yale University\Research\Data\mooring_data\mooring_' + string(mooring_no) + '\uls' + string(yrs{ii,1}) + string(mooring_no) + '_daily*.mat');
load(files.name);

%%% Formatting time
clear tmp_time
for u = 1:size(dates, 1)
    tmp_time(u) = datetime(str2double(dates(u, 1:4)), str2double(dates(u, 6:7)), str2double(dates(u, 9:10)), 'TimeZone', 'UTC');
    if ii == 18
        tmp_time(u) = datetime(dates(u,:), 'InputFormat', 'dd-MMM-uuuu', 'TimeZone', 'UTC');
    end
end
tmp_time = tmp_time(:);

%%% Saving data to structure
uls_tmp.datetime = tmp_time;
uls_tmp.ssh_uncorrected = WL;
ssh_corrected = uls_tmp.ssh_uncorrected;

%%%
figure()
Handle = plot(uls_tmp.datetime, uls_tmp.ssh_uncorrected);
brush on

%%

xd = get(Handle, 'XData');
brush = get(Handle, 'BrushData');
eddy_dates1 = xd(logical(brush));
close all

%%% Removing obvious eddies
[~, ind] = intersect(tmp_time, eddy_dates1);
ssh_corrected(ind) = NaN;

%%% Removing mooring stretch
x = datenum(tmp_time);
y = ssh_corrected;
x = x(:); y = y(:);
x0 = mean(x,'omitnan');      % or median(x,'omitnan')
t  = x - x0;                 % centered days
p = robustfit([t t.^2], y);  % intercept is added automatically
yhat = p(1) + p(2)*t + p(3)*t.^2;
[min_val, ind] = min(yhat);
yhat(ind:end) = min_val;
yd = y - yhat;

%%% Plotting uls data with obvious eddies and stretch removed
figure()
Handle2 = plot(uls_tmp.datetime, yd);
brush on

%%

xd = get(Handle2, 'XData');
brush = get(Handle2, 'BrushData');
eddy_dates2 = xd(logical(brush));

eddy_dates = [eddy_dates1 eddy_dates2];

%%% Saving eddy locations
save('C:\Users\jak279\OneDrive - Yale University\Research\Data\mooring_data\mooring_' + string(mooring_no) + '\uls' + string(yrs{ii,1}) + string(mooring_no) + '_eddy_locations.mat', 'eddy_dates', 'mooring_no')

close all

%%
uu = 1;
for ii = 1:length(yrs)

    %%% Skipping years without data
    if strcmp(mooring_no, 'b') & strcmp(string(yrs{ii,1}), '09')
        uu = uu + 1;
        continue
    elseif strcmp(mooring_no, 'd') & strcmp(string(yrs{ii,1}), '03')
        uu = uu + 1;
        continue
    elseif strcmp(mooring_no, 'd') & strcmp(string(yrs{ii,1}), '04')
        uu = uu + 1;
        continue
    elseif strcmp(mooring_no, 'd') & strcmp(string(yrs{ii,1}), '05')
        uu = uu + 1;
        continue
    end

    %%% Loading data
    files = dir('C:\Users\jak279\OneDrive - Yale University\Research\Data\mooring_data\mooring_' + string(mooring_no) + '\uls' + string(yrs{ii,1}) + string(mooring_no) + '_daily*.mat');
    load(files.name);

    %%% Formatting time
    clear tmp_time
    for u = 1:size(dates, 1)
        tmp_time(u) = datetime(str2double(dates(u, 1:4)), str2double(dates(u, 6:7)), str2double(dates(u, 9:10)), 'TimeZone', 'UTC');
        if ii == 18
            tmp_time(u) = datetime(dates(u,:), 'InputFormat', 'dd-MMM-uuuu', 'TimeZone', 'UTC');
        end
    end
    tmp_time = tmp_time(:);

    %%% Saving data to structure
    uls_tmp.datetime = tmp_time;
    uls_tmp.draft = IDS(:, 2);
    uls_tmp.draft_minimum = IDS(:, 4);
    uls_tmp.draft_maximum = IDS(:, 5);
    uls_tmp.draft_median = IDS(:, 6);
    uls_tmp.ssh_uncorrected = WL;
    ssh_corrected = uls_tmp.ssh_uncorrected;

    %%% Getting dates with eddies, if applicable
    file_name = 'C:\Users\jak279\OneDrive - Yale University\Research\Data\mooring_data\mooring_' + string(mooring_no) + '\uls' + string(yrs{ii,1}) + string(mooring_no) + '_eddy_locations.mat';
    if isfile(file_name)
        load(file_name);
        [~, ind] = intersect(uls_tmp.datetime, eddy_dates);
        ssh_corrected(ind) = NaN;
    end

    %%% Removing mooring stretch
    x = datenum(uls_tmp.datetime);
    y = ssh_corrected;
    x = x(:); y = y(:);
    x0 = mean(x,'omitnan');      % or median(x,'omitnan')
    t  = x - x0;                 % centered days
    p = robustfit([t t.^2], y);  % intercept is added automatically
    yhat = p(1) + p(2)*t + p(3)*t.^2;
    [min_val, ind] = min(yhat);
    yhat(ind:end) = min_val;
    yd = y - yhat;

    %%% Adding SLP data to obtain p_ULS
    ind = slp.datetime >= uls_tmp.datetime(1) & slp.datetime <= uls_tmp.datetime(end);
    slp_tmp = slp.(string(mooring_no))(ind);
    pressure = yd + slp_tmp;

    figure('Position', [10 10 1000 500])
    tiledlayout(1, 2)

    nexttile()
    plot(x, y)
    hold on
    yhat(yhat == 0) = NaN;
    plot(x, yhat)
    ylabel('[m]')

    %%% Computing ssh anomalies
    uls_tmp.ssh_anomaly = yd - mean(yd, 'omitnan');
    uls_tmp.pressure_anomaly = pressure - mean(pressure, 'omitnan');
    
    nexttile()
    plot(uls_tmp.datetime, pressure)
    hold on
    plot(uls_tmp.datetime, uls_tmp.pressure_anomaly)

    %%% Saving start/end dates
    uls_tmp.start = uls_tmp.datetime(1);
    uls_tmp.end = uls_tmp.datetime(end);

    %%% Dividing the 2018-2021 deployment into ~1 year segments
    if strcmp(yrs{ii,1}, '18')
        n = floor(length(tmp_time) / 3);
        dividing_indices = [0, n, 2*n, length(tmp_time)];

        for m = 1:3
            idx = (dividing_indices(m)+1):dividing_indices(m+1);
            uls_tmp.datetime = tmp_time(idx);
            uls_tmp.draft = IDS(idx, 2);
            uls_tmp.draft_minimum = IDS(idx, 4);
            uls_tmp.draft_maximum = IDS(idx, 5);
            uls_tmp.draft_median = IDS(idx, 6);

            uls_tmp.ssh_uncorrected = WL(idx);
            uls_tmp.ssh_anomaly = yd(idx) - mean(yd(idx), 'omitnan');
            uls_tmp.pressure_anomaly = pressure(idx) - mean(pressure(idx), 'omitnan');
            uls_tmp.start = uls_tmp.datetime(1);
            uls_tmp.end = uls_tmp.datetime(end);

            uls_all(uu,1) = uls_tmp;
            uu = uu + 1;
        end
        clear uls_tmp idx n m dividing_indices ssh_corrected mooring_stretch p s ins ssh_no_eddies IDS IS BETA BTBETA dates ii WL WLS yday OWBETA ind u tmp_time
        continue
    end

    %%% Saving structure
    uls_all(uu,1) = uls_tmp;
    uu = uu + 1;
    clear uls_tmp idx n m dividing_indices ssh_corrected mooring_stretch p s ins ssh_no_eddies IDS IS BETA BTBETA dates ii WL WLS yday OWBETA ind u tmp_time
end

%%


%%% Combining data into one long time series
clear combined
uls = uls_all;
combined.datetime = uls(1).datetime;
combined.draft = uls(1).draft;
combined.draft_minimum = uls(1).draft_minimum;
combined.draft_maximum = uls(1).draft_maximum;
combined.draft_median = uls(1).draft_median;
combined.ssh_uncorrected = uls(1).ssh_uncorrected;
combined.ssh_anomaly = uls(1).ssh_anomaly;
combined.pressure_anomaly = uls(1).pressure_anomaly;
if size(uls,1) > 1
    for i = 2:size(uls, 1)
        combined.datetime = [combined.datetime; uls(i).datetime];
        combined.draft = [combined.draft; uls(i).draft];
        combined.draft_minimum = [combined.draft_minimum; uls(i).draft_minimum];
        combined.draft_maximum = [combined.draft_maximum; uls(i).draft_maximum];
        combined.draft_median = [combined.draft_median; uls(i).draft_median];
        combined.ssh_uncorrected = [combined.ssh_uncorrected; uls(i).ssh_uncorrected];
        combined.ssh_anomaly = [combined.ssh_anomaly; uls(i).ssh_anomaly];
        combined.pressure_anomaly = [combined.pressure_anomaly; uls(i).pressure_anomaly];
    end
end
clear i

%%% Add NaNs to fill data gaps
combined_final.datetime = tgrid;
combined_final.draft = NaN(size(tgrid));
combined_final.draft_minimum = NaN(size(tgrid));
combined_final.draft_maximum = NaN(size(tgrid));
combined_final.draft_median = NaN(size(tgrid));
combined_final.ssh_uncorrected = NaN(size(tgrid));
combined_final.ssh_anomaly = NaN(size(tgrid));
combined_final.pressure_anomaly = NaN(size(tgrid));
for i = 1:length(combined_final.datetime)
    [~,idx] = ismember(combined_final.datetime(i), combined.datetime);
    if idx == 0
        continue
    end
    combined_final.draft(i) = combined.draft(idx);
    combined_final.draft_minimum(i) = combined.draft_minimum(idx);
    combined_final.draft_maximum(i) = combined.draft_maximum(idx);
    combined_final.draft_median(i) = combined.draft_median(idx);
    combined_final.ssh_uncorrected(i) = combined.ssh_uncorrected(idx);
    combined_final.ssh_anomaly(i) = combined.ssh_anomaly(idx);
    combined_final.pressure_anomaly(i) = combined.pressure_anomaly(idx);
end

%%% Adding column with start/end dates
combined_final.start = NaT(1, length(uls), 'TimeZone', 'UTC');
combined_final.end = NaT(1, length(uls), 'TimeZone', 'UTC');
for i = 1:length(uls)
    if isempty(uls(i).start)
        continue
    end
    combined_final.start(i) = uls(i).start;
    combined_final.end(i) = uls(i).end;
end

uls = combined_final;
uls.lat = lat;
uls.lon = lon;
clear combined_final idx


%%% Saving data
dt_string = char(dt);
dt_string = string(dt_string(~isspace(dt_string)));
save('C:\Users\jak279\OneDrive - Yale University\Research\Data\OBP_interp\ULS_mooring_' + string(mooring_no) + '_' + string(dt_string) + '.mat',  'uls', 'uls_all', '-v7.3');
clear dt_string

%%

figure()
hold on
yyaxis right
%plot(uls.datetime, uls.pressure)
yyaxis left
plot(uls.datetime, uls.pressure_anomaly)
