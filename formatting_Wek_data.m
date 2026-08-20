clc; clear;

%%% Mooring number
mooring_no = 'a';
if strcmp(mooring_no, 'a')
    target_lat = 75;
    target_lon = -150;
elseif strcmp(mooring_no, 'b')
    target_lat = 78;
    target_lon = -150;
elseif strcmp(mooring_no, 'd')
    target_lat = 74;
    target_lon = -140;
end

%%% Loading data
time = ncread('C:\Users\jak279\OneDrive - Yale University\Research\Data\ArcticEkmanPumping.nc', 'time');
datetime_grd = datetime(datenum(2003, 1, 5) + time, 'ConvertFrom', 'datenum', 'TimeZone', 'UTC');
lat_grd = ncread('C:\Users\jak279\OneDrive - Yale University\Research\Data\ArcticEkmanPumping.nc', 'lat');
lon_grd = ncread('C:\Users\jak279\OneDrive - Yale University\Research\Data\ArcticEkmanPumping.nc', 'lon');
we = ncread('C:\Users\jak279\OneDrive - Yale University\Research\Data\ArcticEkmanPumping.nc', 'we');
alphaweice = ncread('C:\Users\jak279\OneDrive - Yale University\Research\Data\ArcticEkmanPumping.nc', 'alphaweice');
alphaweair = ncread('C:\Users\jak279\OneDrive - Yale University\Research\Data\ArcticEkmanPumping.nc', 'alphaweair');
rho0 = 1027.5;
f0 = 1.46e-4;

%%% Getting distance between moorings and grid points
wgs84 = wgs84Ellipsoid("km");
dist_2_mooring = distance(target_lat, target_lon, lat_grd, lon_grd, wgs84);
radius_for_averaging = 100;
indk = find(dist_2_mooring < radius_for_averaging);

%%% Taking stresses closest to each mooring
wek.datetime = datetime_grd;

wek.target_lat = target_lat;
wek.target_lon = target_lon;
wek.radius_km = radius_for_averaging;
wek.we = mean(we(indk,:), 1, 'omitnan');
wek.we = wek.we(:);
wek.alphaweice = mean(alphaweice(indk,:));
wek.alphaweice = wek.alphaweice(:);
wek.alphaweair = mean(alphaweair(indk,:));
wek.alphaweair = wek.alphaweair(:);

%%% Getting area of grid cell
cell_area = (25e3)^2; % 25km Equal-Area Scalable Earth (EASE) grid
wek.mass_flux = wek.we * cell_area;

%%% Saving data
dt = days(7);
dt_string = char(dt);
dt_string = string(dt_string(~isspace(dt_string)));
save('C:\Users\jak279\OneDrive - Yale University\Research\Data\OBP_interp\Wek_mooring_' + string(mooring_no) + '_' + string(dt_string) + '.mat', 'wek')

%%% Plotting
figure()
hold on
plot(wek.datetime, wek.we, 'DisplayName', 'Ekman Pumping')
%plot(wek.datetime, wek.alphaweice, 'DisplayName', 'Ice-Ocean Ekman Pumping')
%plot(wek.datetime, wek.alphaweair, 'DisplayName', 'Wind-Ocean Ekman Pumping')
legend()
grid on
ylabel('Ekman Pumping [m/s]')

%% Map figure showing the area averaged over

land = readgeotable("landareas.shp");
landColor = [0.8 0.8 0.8];

%%% Creating figure
figure('Position', [10 10 700 700]);
set(gcf, 'Color', 'w')

%%% Plotting Ekman pumping at a single time slice
ax_main = axesm('polycon', 'MapLatLimit', [70.5 80.5], 'MapLonLimit', [-170 -130], 'Frame', 'on', 'Grid', 'on', 'MeridianLabel', 'on', 'ParallelLabel', 'on', 'PLineLocation', 5, 'MLineLocation', 20);
scatterm(lat_grd, lon_grd, 35, we(:,2), 'filled', 'Marker', 'square')
hold on
scatterm(lat_grd(indk), lon_grd(indk), 35, we(indk,2), 'filled', 'MarkerFaceColor', 'r')
scatterm(target_lat, target_lon, 50, 'filled', 'MarkerFaceColor', 'k')
set(ax_main, 'Color', 'none', 'Box', 'off', 'XColor', 'none', 'YColor', 'none');
colormap(cmocean('balance')); 
cb = colorbar;
clim([-5e-6 5e-6])
cb.Position = [0.89, 0.2, 0.03, 0.6];
geoshow(ax_main, land, FaceColor = landColor);


