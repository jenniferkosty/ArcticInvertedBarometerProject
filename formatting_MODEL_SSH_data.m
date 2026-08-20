clc; clear;

%% Formatting model output data

%%% Mooring location
mooring_no = 'np';

if strcmp(mooring_no, 'a')
    lat = 75;
    lon = -150;
    col = 6;
elseif strcmp(mooring_no, 'b')
    lat = 78;
    lon = -150;
    col = 7;
elseif strcmp(mooring_no, 'd')
    lat = 74;
    lon = -140;
    col = 9;
elseif strcmp(mooring_no, 'np')
    lat = 89.25;
    lon = 60;
    col = 5;
end

%%% Loading data
u = 3;
labels_all = {'Control', 'SLP', 'Wind'};
label = labels_all{u};

if strcmp(label, 'Control')
    file_label = 'CTR';
elseif strcmp(label, 'SLP')
    file_label = 'SLP';
elseif strcmp(label, 'Wind')
    file_label = 'WIND';
end

%%% Moorings
file_dir = 'C:\Users\jak279\OneDrive - Yale University\Research\Data\ModelOutput081126\';
T = readtable(string(file_dir) + 'HOURLY-ssh-moorings-ice-' + string(file_label) + '-Aug-12');

year1 = 2003;
year2 = 2023;

%%% Formatting
ssh_tmp.ssh = T.('Var' + string(col));
year = T.Var1;
month = T.Var2;
day = T.Var3;
hour = T.Var4;
ssh_tmp.datetime = datetime(year, month, day, hour, 0, 0, 'TimeZone', 'UTC');
clear year month day hour T i

%%% Dividing into year-long segments
first_year = year1;
last_year = year2;
dividing_dates = datetime(first_year:last_year+1, 1, 1, 'TimeZone', 'UTC');
clear ssh_highres_all
for i = 1:length(dividing_dates)-1
    ind = ssh_tmp.datetime >= dividing_dates(i) & ssh_tmp.datetime < dividing_dates(i+1);
    ssh_highres_all(i,1).datetime = ssh_tmp.datetime(ind);
    ssh_highres_all(i,1).ssh = ssh_tmp.ssh(ind);
    ssh_highres_all(i,1).ssh_anomaly = ssh_highres_all(i).ssh - mean(ssh_highres_all(i).ssh, 'omitnan');
    ssh_highres_all(i,1).start = ssh_highres_all(i).datetime(1);
    ssh_highres_all(i,1).end = ssh_highres_all(i).datetime(end);
end
clear dividing_dates i

% %%% Detiding annual segments using t-tide
% for i = 1:length(ssh_highres_all)
%     if ~isempty(ssh_highres_all(i).ssh_anomaly)
% 
%         t = ssh_highres_all(i).datetime;
%         y = ssh_highres_all(i).ssh_anomaly;
%         t_num = datenum(t);
%         ttide_dt = (t_num(2) - t_num(1)) * 24; % units of hours
% 
%         %%% Finding NaNs in segment
%         nan_mask = ~isfinite(y);
%         gap_starts = find(diff([0; nan_mask]) == 1);
%         gap_ends   = find(diff([nan_mask; 0]) == -1);
%         gap_lengths_hours = (gap_ends - gap_starts + 1) * ttide_dt;
% 
%         %%% Skipping too many NaNs
%         if sum(gap_lengths_hours) > 0.5 * ((t_num(end) - t_num(1)) * 24)
%             detided = NaN(size(y));
%             continue
%         end
% 
% 
%         %%% Detiding
%         [~, ~, ~, pout] = t_tide(y, ...
%             'interval', ttide_dt, ...
%             'start',    t_num(1), ...
%             'latitude', lat);
%         detided = y - pout;
% 
%         %%% Ensure long gap positions stay NaN 
%         detided(nan_mask) = NaN;
% 
%         %%% Saving
%         ssh_highres_all(i).detided = detided;
% 
%     end
% 
%     clear detided good seg_starts seg_ends short_gaps idx y_filled g nan_mask gap_starts gap_ends gap_length_hours t_num
% end

%%% Creating time grid
dt = hours(1);
tgrid = datetime(first_year, 1, 1, "TimeZone", "UTC"):dt:datetime(last_year+1, 1, 1, "TimeZone", "UTC");
tgrid = tgrid(1:end-1);
tgrid = tgrid(:);

%%% Combining high-res data into one long time series
clear combined
ssh = ssh_highres_all;
combined.datetime = ssh(1).datetime;
combined.ssh = ssh(1).ssh;
combined.ssh_anomaly = ssh(1).ssh_anomaly;
%combined.detided = ssh(i).detided;
if size(ssh_highres_all,1) > 1
    for i = 2:size(ssh_highres_all, 1)
        combined.datetime = [combined.datetime; ssh(i).datetime];
        combined.ssh = [combined.ssh; ssh(i).ssh];
        combined.ssh_anomaly = [combined.ssh_anomaly; ssh(i).ssh_anomaly];
        %combined.detided = [combined.detided; ssh(i).detided];
    end
end
clear i

%%% Add NaNs to fill data gaps
combined_final.datetime = tgrid;
combined_final.ssh = NaN(size(tgrid));
combined_final.ssh_anomaly = NaN(size(tgrid));
%combined_final.detided = NaN(size(tgrid));
for i = 1:length(combined_final.datetime)
    [~,idx] = ismember(combined_final.datetime(i), combined.datetime);
    if idx == 0
        continue
    end
    combined_final.ssh(i) = combined.ssh(idx);
    combined_final.ssh_anomaly(i) = combined.ssh_anomaly(idx);
    %combined_final.detided(i) = combined.detided(idx);
end

%%% Adding column with start/end dates
combined_final.start = NaT(1, length(ssh), 'TimeZone', 'UTC');
combined_final.end = NaT(1, length(ssh), 'TimeZone', 'UTC');
for i = 1:length(ssh)
    if isempty(ssh(i).start)
        continue
    end
    combined_final.start(i) = ssh(i).start;
    combined_final.end(i) = ssh(i).end;
end

ssh_highres = combined_final;
ssh_highres.lat = lat;
ssh_highres.lon = lon;
clear combined_final idx

%%% Saving data
dt_string = char(dt);
dt_string = string(dt_string(~isspace(dt_string)));
ssh = ssh_highres;
ssh_all = ssh_highres_all;
save('C:\Users\jak279\OneDrive - Yale University\Research\Data\OBP_interp\ModelOutput' + string(label) + '_mooring_' + string(mooring_no) + '_' + string(dt_string) + '.mat', 'ssh', 'ssh_all')
clear dt_string ssh ssh_all 


%%

%%% Creating coarser time grid
dt = hours(6);
tgrid = datetime(first_year, 1, 1, "TimeZone", "UTC"):dt:datetime(last_year+1, 1, 1, "TimeZone", "UTC");
tgrid = tgrid(1:end-1);
tgrid = tgrid(:);

%%% Running mean 
clear ssh_all
for ii = 1:size(ssh_highres_all,1)

    %%% Getting time array for each year
    t1 = datetime(year(ssh_highres_all(ii).start), 1, 1, "TimeZone", "UTC");
    t2 = datetime(year(ssh_highres_all(ii).end)+1, 12, 31, "TimeZone", "UTC");
    t_arr.datetime = tgrid(tgrid >= t1 & tgrid <= t2);

    %%% Selecting data for 1 year
    ssh = ssh_highres_all(ii,1);

    %%% Getting data on time grid
    clear ssh_interp
    ssh_interp.ssh = NaN(length(t_arr.datetime),1);
    ssh_interp.ssh_anomaly = NaN(length(t_arr.datetime),1);
    for i = 2:length(t_arr.datetime) - 1
        ind = (datenum(ssh.datetime) > (datenum(t_arr.datetime(i)) - (dt/2))) & (datenum(ssh.datetime) <= (datenum(t_arr.datetime(i)) + (dt/2)));
        ssh_interp.ssh(i) = mean(ssh.ssh(ind));
        ssh_interp.ssh_anomaly(i) = mean(ssh.ssh_anomaly(ind));
    end
    ssh_interp.datetime = t_arr.datetime;
    ssh_interp.start = ssh.start;
    ssh_interp.end = ssh.end;

    %%% Removing NaNs
    idx = ~isnan(ssh_interp.ssh);
    ssh_interp.ssh = ssh_interp.ssh(idx);
    ssh_interp.ssh_anomaly = ssh_interp.ssh_anomaly(idx);
    ssh_interp.datetime = ssh_interp.datetime(idx);
    ssh_all(ii,1) = ssh_interp;
end

%%% Combining daily data into one long time series
clear combined
ssh = ssh_all;
combined.ssh = ssh(1).ssh;
combined.datetime = ssh(1).datetime;
combined.ssh_anomaly = ssh(1).ssh_anomaly;
if size(ssh_all,1) > 1
    for i = 2:size(ssh_all, 1)
        combined.ssh = [combined.ssh; ssh(i).ssh];
        combined.ssh_anomaly = [combined.ssh_anomaly; ssh(i).ssh_anomaly];
        combined.datetime = [combined.datetime; ssh(i).datetime];
    end
end
clear i

%%% Add NaNs to fill data gaps
combined_final.datetime = tgrid;
combined_final.ssh = NaN(size(tgrid));
combined_final.ssh_anomaly = NaN(size(tgrid));
for i = 1:length(combined_final.datetime)
    [~,idx] = ismember(combined_final.datetime(i), combined.datetime);
    if idx == 0
        continue
    end
    combined_final.ssh(i) = combined.ssh(idx);
    combined_final.ssh_anomaly(i) = combined.ssh_anomaly(idx);
end

%%% Adding column with start/end dates
combined_final.start = NaT(1, length(ssh), 'TimeZone', 'UTC');
combined_final.end = NaT(1, length(ssh), 'TimeZone', 'UTC');
for i = 1:length(ssh)
    if isempty(ssh(i).start)
        continue
    end
    combined_final.start(i) = ssh(i).start;
    combined_final.end(i) = ssh(i).end;
end

clear ii i t1 t2 t_arr ind

%%% Saving data
ssh = combined_final;
ssh.lat = lat;
ssh.lon = lon;
dt_string = char(dt);
dt_string = string(dt_string(~isspace(dt_string)));
save('C:\Users\jak279\OneDrive - Yale University\Research\Data\OBP_interp\ModelOutput' + string(label) + '_mooring_' + string(mooring_no) + '_' + string(dt_string) + '.mat', 'ssh', 'ssh_all')






