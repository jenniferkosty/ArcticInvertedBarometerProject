%% Script to produce wind-stress curl (WSC) time series for a specified lat/lon
clear; clc;

%%% Setting location
mooring_no = 'np';

if strcmp(mooring_no, 'a')
    target_lat = 75;
    target_lon = -150;
elseif strcmp(mooring_no, 'b')
    target_lat = 78;
    target_lon = -150;
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

%%% Creating time grid (data have daily resolution)
dt = days(1);
first_year = 2003;
last_year = 2025;
tgrid = datetime(first_year, 1, 1, "TimeZone", "UTC"):dt:datetime(last_year+1, 1, 1, "TimeZone", "UTC");
tgrid = tgrid(1:end-1);
tgrid = tgrid(:);

%%% Getting years of time array
yrs_all = unique(year(tgrid));

%%% Formatting u-wind data
for ii = 1:length(yrs_all)
    nc_name = 'C:\Users\jak279\OneDrive - Yale University\Research\Data\2D Data Products\ERA5 Winds\era5_u_' + string(yrs_all(ii)) + '.nc';
    data_all{ii} = ncread(nc_name, 'u10'); 
    time_all{ii} = datetime(ncread(nc_name, 'valid_time') + datenum(yrs_all(ii), 1, 1), 'ConvertFrom', 'datenum', 'TimeZone', 'UTC'); 
end
full.data = data_all{1};
full.datetime = time_all{1};
if length(data_all) > 1
    for i = 2:length(data_all)
        full.data = cat(3, full.data, data_all{i});
        full.datetime = cat(1, full.datetime, time_all{i});
    end
end
lati = ncread(nc_name, 'latitude');
loni = ncread(nc_name, 'longitude'); 
ind = loni < 0;
loni(ind) = loni + 360;
[full.lat, full.lon] = meshgrid(lati,loni); % transforms vectors lati and loni into arrays 
u = full.data;

%%% Finding grid cell closest to chosen lat/lon
isLat = full.lat > lat-3 & full.lat < lat+3; % find latitude indices around target point
isLon = full.lon > lon-3 & full.lon < lon+3; % find longitude indices around target point
nearI = find(isLat & isLon);
nNear = length(nearI); % all the indices in the vicinity
dist = distance(full.lat(nearI),full.lon(nearI),repmat(lat,nNear,1),repmat(lon,nNear,1)); % compute distance from all the indices in the vicinity
[~,I] = min(dist);
closestI = nearI(I);
[r,c] = ind2sub(size(full.lat), closestI); % index of the closest grid cell

%%% Creating wind time series (daily resolution)
wind_highres.u_speed = squeeze(u(r,c,:));
wind_highres.datetime = full.datetime;
wind_highres.target_lat = target_lat;
wind_highres.target_lon = target_lon;
wind_highres.actual_lat = full.lat(r,c);
wind_highres.actual_lon = full.lon(r,c);
if wind_highres.actual_lon > 180
    wind_highres.actual_lon = wind_highres.actual_lon - 360;
end

%%% Formatting v-wind data
for ii = 1:length(yrs_all)
    nc_name = 'C:\Users\jak279\OneDrive - Yale University\Research\Data\2D Data Products\ERA5 Winds\era5_v_' + string(yrs_all(ii)) + '.nc';
    data_all{ii} = ncread(nc_name, 'v10'); 
    time_all{ii} = datetime(ncread(nc_name, 'valid_time') + datenum(yrs_all(ii), 1, 1), 'ConvertFrom', 'datenum', 'TimeZone', 'UTC'); 
end
full.data = data_all{1};
full.datetime = time_all{1};
if length(data_all) > 1
    for i = 2:length(data_all)
        full.data = cat(3, full.data, data_all{i});
        full.datetime = cat(1, full.datetime, time_all{i});
    end
end
lati = ncread(nc_name, 'latitude');
loni = ncread(nc_name, 'longitude'); % degrees east: [-180:180]
ind = loni < 0;
loni(ind) = loni + 360;
[full.lat, full.lon] = meshgrid(lati,loni); % transforms vectors lati and loni into arrays 
v = full.data;

%%% Finding grid cell closest to chosen lat/lon
isLat = full.lat > lat-3 & full.lat < lat+3; % find latitude indices around target point
isLon = full.lon > lon-3 & full.lon < lon+3; % find longitude indices around target point
nearI = find(isLat & isLon);
nNear = length(nearI); % all the indices in the vicinity
dist = distance(full.lat(nearI),full.lon(nearI),repmat(lat,nNear,1),repmat(lon,nNear,1)); % compute distance from all the indices in the vicinity
[~,I] = min(dist);
closestI = nearI(I);
[r,c] = ind2sub(size(full.lat), closestI); % index of the closest grid cell

%%% Creating wind time series (daily resolution)
wind_highres.v_speed = squeeze(v(r,c,:));
wind_highres.datetime = full.datetime;
wind_highres.target_lat = target_lat;
wind_highres.target_lon = target_lon;
wind_highres.actual_lat = full.lat(r,c);
wind_highres.actual_lon = full.lon(r,c);
if wind_highres.actual_lon > 180
    wind_highres.actual_lon = wind_highres.actual_lon - 360;
end

%%% Computing curl
clear full
[full.lat, full.lon] = meshgrid(lati,loni);
full.datetime = tgrid;
rho0 = 1.25;
cd = 0.00125;
del = 1;
taux = rho0 .* cd .* u.^2;
tauy = rho0 .* cd .* v.^2;
curl_arr = NaN(size(u));
for i =  2:length(loni)-1
    for j = 2:length(lati)-1
        dtauy = tauy(i+del,j,:) - tauy(i-del,j,:);
        x_m = gc_dist(full.lat(i+del, j), full.lon(i+del,j), full.lat(i-del,j), full.lon(i-del,j)) * 1000;
        dtauy_dx = dtauy ./ x_m;

        dtaux = taux(i,j+del,:) - taux(i,j-del,:);
        y_m = gc_dist(full.lat(i,j+del), full.lon(i,j+del), full.lat(i,j-del), full.lon(i,j-del)) * 1000;
        dtaux_dy = dtaux ./ y_m;

        dtauy_dx = reshape(dtauy_dx,[],1);
        dtaux_dy = reshape(dtaux_dy,[],1);

        curl_arr(i,j,:) = dtauy_dx - dtaux_dy;
    end
end
wind_stress_curl = -curl_arr;

%%% Finding grid cell closest to chosen lat/lon
isLat = full.lat > lat-3 & full.lat < lat+3; % find latitude indices around target point
isLon = full.lon > lon-3 & full.lon < lon+3; % find longitude indices around target point
nearI = find(isLat & isLon);
nNear = length(nearI); % all the indices in the vicinity
e = referenceEllipsoid('WGS84', 'km');
dist = distance(full.lat(nearI),full.lon(nearI),repmat(lat,nNear,1),repmat(lon,nNear,1), e); % compute distance from all the indices in the vicinity
%[~,I] = min(dist);

%%% Averaging over values near target lat/lon
radius_for_averaging = 100;
I = dist < radius_for_averaging;
closestI = nearI(I);
[r,c] = ind2sub(size(full.lat), closestI); % index of the closest grid cell

%%% Creating wind time series (daily resolution)
wind_highres.curl = squeeze(mean(wind_stress_curl(r,c,:), [1, 2]));

%%% Using running mean to obtain desired time grid
wind = wind_highres;
if dt > days(1)
    wind.datetime = tgrid;
    wind.u_speed = NaN(size(tgrid));
    wind.v_speed = NaN(size(tgrid));
    wind.curl = NaN(size(tgrid));
    % for i = 1:length(tgrid)
    %     ind = (wind_highres.datetime > tgrid(i) - (dt/2)) & (wind_highres.datetime <= tgrid(i) + (dt/2));
    %     wind.u(i) = mean(wind_highres.u(ind));
    %     wind.v(i) = mean(wind_highres.v(ind));
    %     wind.curl(i) = mean(wind_highres.curl(ind));
    % end

    [wind.datetime, wind.u_speed] = running_mean(wind_highres.datetime, wind_highres.u_speed, tgrid);
    [wind.datetime, wind.v_speed] = running_mean(wind_highres.datetime, wind_highres.v_speed, tgrid);
    [wind.datetime, wind.curl] = running_mean(wind_highres.datetime, wind_highres.curl, tgrid);
end

%%% Saving data
dt_string = char(dt);
dt_string = string(dt_string(~isspace(dt_string)));
save('C:\Users\jak279\OneDrive - Yale University\Research\Data\OBP_interp\WSC_mooring_' + string(mooring_no) + '_' + string(dt_string) + '.mat', 'wind');


