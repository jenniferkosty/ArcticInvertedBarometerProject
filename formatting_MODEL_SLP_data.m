clc; clear;

%%% Mooring location
mooring_no = 'np';

if strcmp(mooring_no, 'a')
    lat = 75;
    lon = -150;
elseif strcmp(mooring_no, 'b')
    lat = 78;
    lon = -150;
elseif strcmp(mooring_no, 'd')
    lat = 74;
    lon = -140;
elseif strcmp(mooring_no, 'np')
    lat = 89.25;
    lon = 60;
end

%%% Loading data
T = readtable('C:\Users\jak279\OneDrive - Yale University\Research\Data\AndreyModelOutput-2003-2023\SLP-6-hourly-2003-2023-AT-MOORING-and-NP');

%%% Formatting
slp_highres.pressure = T.(string(upper(mooring_no))) ./ 100; % converting millibars to dbars
year = T.Year;
month = T.Month;
day = T.Day;
hour = T.hour;
slp_highres.datetime = datetime(year, month, day, hour, 0, 0, 'TimeZone', 'UTC');
clear year month day hour T i

%%% Creating time grid
dt = hours(6);
first_year = 2003;
last_year = 2023;
tgrid = datetime(first_year, 1, 1, "TimeZone", "UTC"):dt:datetime(last_year+1, 1, 1, "TimeZone", "UTC");
tgrid = tgrid(1:end-1);
tgrid = tgrid(:);

%%% Using running mean to obtain desired time grid
slp = slp_highres;
if dt > hours(1)
    slp.datetime = tgrid;
    slp.pressure = NaN(size(tgrid));
    for i = 1:length(tgrid)
        ind = (slp_highres.datetime > tgrid(i) - (dt/2)) & (slp_highres.datetime <= tgrid(i) + (dt/2));
        slp.pressure(i) = mean(slp_highres.pressure(ind));
    end
end

%%% Saving data
dt_string = char(dt);
dt_string = string(dt_string(~isspace(dt_string)));
save('C:\Users\jak279\OneDrive - Yale University\Research\Data\OBP_interp\ModelSLP_mooring_' + string(mooring_no) + '_' + string(dt_string) + '.mat', 'slp');
clear dt_string

