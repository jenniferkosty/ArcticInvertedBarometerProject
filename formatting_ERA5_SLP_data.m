%% Script to produce SLP time series for a specified lat/lon
clear; clc;

%%% Setting location
mooring_no = 'np';

if strcmp(mooring_no, 'a')
    target_lat = 75;
    target_lon = -150;
elseif strcmp(mooring_no, 'b')
    target_lat = 78;
    target_lon = -150;
elseif strcmp(mooring_no, 'c')
    target_lat = 77;
    target_lon = -140;
elseif strcmp(mooring_no, 'd')
    target_lat = 74;
    target_lon = -140;
elseif strcmp(mooring_no, 'np')
    target_lat = 89.25;
    target_lon = 60.35;
end

%%% Reformatting lon value 
lat = target_lat;
lon = target_lon;
if lon < 0
    lon = lon + 360;
end

%%% Creating time grid
dt = hours(1);
first_year = 2024;
last_year = 2024;
tgrid = datetime(first_year, 1, 1, "TimeZone", "UTC"):dt:datetime(last_year+1, 1, 1, "TimeZone", "UTC");
tgrid = tgrid(1:end-1);
tgrid = tgrid(:);

%%% Getting years of time array
yrs_all = unique(year(tgrid));

%%% Loading SLP data
for ii = 1:length(yrs_all)
    nc_name = 'C:\Users\jak279\OneDrive - Yale University\Research\Data\2D Data Products\ERA5 SLP\era5_slp_' + string(yrs_all(ii)) + '.nc';
    slpdata_all{ii} = double(ncread(nc_name, 'msl')./10000); % /10000 to convert Pascals to dbar
    slptime_all{ii} = datetime(ncread(nc_name, 'valid_time'), 'ConvertFrom', 'posixtime', 'TimeZone', 'UTC'); % seconds since 1970-01-01
end
slp_full.data = slpdata_all{1};
slp_full.datetime = slptime_all{1};
if length(slpdata_all) > 1
    for i = 2:length(slpdata_all)
        slp_full.data = cat(3, slp_full.data, slpdata_all{i});
        slp_full.datetime = cat(1, slp_full.datetime, slptime_all{i});
    end
end
lati = double(ncread(nc_name, 'latitude'));
loni = double(ncread(nc_name, 'longitude')); 
ind = loni < 0;
loni(ind) = loni(ind) + 360;
[slp_full.lat, slp_full.lon] = meshgrid(lati,loni); % grid of lat/lon
clear lati loni ii

%%% Finding grid cell closest to chosen lat/lon
isLat = slp_full.lat > lat-3 & slp_full.lat < lat+3; % find latitude indices around target point
isLon = slp_full.lon > lon-3 & slp_full.lon < lon+3; % find longitude indices around target point
nearI = find(isLat & isLon);
nNear = length(nearI); % all the indices in the vicinity
dist = distance(slp_full.lat(nearI),slp_full.lon(nearI),repmat(lat,nNear,1),repmat(lon,nNear,1)); % compute distance from all the indices in the vicinity
[~,I] = min(dist);
closestI = nearI(I);
[r,c] = ind2sub(size(slp_full.lat), closestI); % index of the closest grid cell

%%% Creating SLP time series (6-hourly resolution)
slp_highres.pressure = squeeze(slp_full.data(r,c,:));
slp_highres.datetime = slp_full.datetime;
slp_highres.target_lat = target_lat;
slp_highres.target_lon = target_lon;
slp_highres.actual_lat = slp_full.lat(r,c);
slp_highres.actual_lon = slp_full.lon(r,c);
if slp_highres.actual_lon > 180
    slp_highres.actual_lon = slp_highres.actual_lon - 360;
end
clear slpdata slpday isLat isLon nearI nNear dist I closestI T loni lati lon lat time slp_full slptime_all i ncname yrs_all

%%% Saving data
slp = slp_highres;
dt = hours(slp_highres.datetime(2) - slp_highres.datetime(1));
dt_string = num2str(dt) + string('hr');
save('C:\Users\jak279\OneDrive - Yale University\Research\Data\OBP_interp\ERA5_SLP_mooring_' + string(mooring_no) + '_' + string(dt_string) + '.mat', 'slp');