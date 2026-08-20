%% Script to produce SIC time series for a specified lat/lon
clear; clc;

%%% Setting location
mooring_no = 'd';

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


%%% Creating time grid (SIC data have daily resolution)
dt = days(1);
first_year = 2003;
last_year = 2025;
tgrid = datetime(first_year, 1, 1, "TimeZone", "UTC"):dt:datetime(last_year+1, 1, 1, "TimeZone", "UTC");
tgrid = tgrid(1:end-1);
tgrid = tgrid(:);

%%% Creating SIC time series
seaicedata = NaN(size(tgrid));
for i = 1:length(tgrid) 
   disp(tgrid(i))

   %%% Extracting time of mooring profile
   mooring_date = datestr(datenum(tgrid(i)), 'yyyymmdd');
   mooring_yr = mooring_date(1:4);
   
   %%% Loading SIC data
   directory = 'C:\Users\jak279\OneDrive - Yale University\Research\Data\2D Data Products\SIC\' + string(mooring_yr) + '\';
   file_info = dir(fullfile(directory, '*' + string(mooring_date)' + '*'));
   nc_name = [file_info.folder '\' file_info.name];

   if str2double(mooring_yr) < 2024
       sea_ice.concentration = ncread(nc_name, 'nsidc_nt_seaice_conc');
       sea_ice.projection = projcrs(3411,"Authority","EPSG");
       sea_ice.x = ncread(nc_name, 'xgrid');
       sea_ice.y = ncread(nc_name, 'ygrid');
   else 
       sea_ice.concentration = ncread(nc_name, 'cdr_seaice_conc');
       sea_ice.projection = projcrs(3411,"Authority","EPSG");
       sea_ice.x = ncread(nc_name, 'x');
       sea_ice.y = ncread(nc_name, 'y');
   end
   [xgrid, ygrid] = meshgrid(sea_ice.x, sea_ice.y);
   [sea_ice.lat, sea_ice.lon] = projinv(sea_ice.projection,double(xgrid), double(ygrid));
   sea_ice.datetime = datetime(mooring_date, 'InputFormat', 'yyyyMMdd', 'TimeZone', 'UTC');

   %%% Finding closest lat/lon for interpolation
   isLat = sea_ice.lat > target_lat-3 & sea_ice.lat < target_lat+3; % find latitude indices around mooring
   isLon = sea_ice.lon > target_lon-3 & sea_ice.lon < target_lon+3; % find longitude indices around mooring
   nearI = find(isLat & isLon);
   nNear = length(nearI); % all the indices in the vicinity
   dist = distance(sea_ice.lat(nearI),sea_ice.lon(nearI),repmat(target_lat,nNear,1),repmat(target_lon,nNear,1)); % compute distance from all the indices in the vicinity
   [~,I] = min(dist);
   closestI = nearI(I);
   [r,c] = ind2sub(size(sea_ice.lat), closestI); %index of the closest grid cell

   % SIC at closest ind
   SIC_data = sea_ice.concentration';
   seaicedata(i) = SIC_data(r,c);
   latdata(i) = sea_ice.lat(r,c);
   londata(i) = sea_ice.lon(r,c);
end
sic.concentration = NaN(size(tgrid));
sic.concentration(1:length(seaicedata)) = seaicedata .*100;
sic.datetime = tgrid;
sic.target_lat = target_lat;
sic.target_lon = target_lon;
sic.actual_lat = latdata;
sic.actual_lon = londata;
clear slpdata slpday c closestI dist l isLat isLon lat lon mooring_date mooring_yr nearI nNear r xgrid ygrid SIC_data seaicedata latdata londata

%%% Saving data
dt_string = char(dt);
dt_string = string(dt_string(~isspace(dt_string)));
save('C:\Users\jak279\OneDrive - Yale University\Research\Data\OBP_interp\SIC_mooring_' + string(mooring_no) + '_' + string(dt_string) + '.mat', 'sic')