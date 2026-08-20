function slp_full = loadingGriddedSLPdata(yrs)

%%% Function loads in the full gridded SLP dataset

% Input: yrs, an array of years, e.g. 2003:2024

% Output: slp_full, a structure with fields for slp, lat, lon, time, etc. 

%%% Loading data for each year
for ii = 1:length(yrs)
    nc_name = "C:\Users\jak279\OneDrive - Yale University\Research\Data\2D Data Products\NCEP SLP\mslp." + string(yrs(ii)) + ".nc";
    slpdata_all{ii} = ncread(nc_name, 'mslp')./10000; % /10000 to convert Pascals to dbar
    slptime_all{ii} = ncread(nc_name, 'time'); % hours since 1800-01-01 00:00:0.0
end
slp_full.pressure = slpdata_all{1};
slp_full.time = slptime_all{1};
if length(slpdata_all) > 1
    for i = 2:length(slpdata_all)
        slp_full.pressure = cat(3, slp_full.pressure, slpdata_all{i});
        slp_full.time = cat(1, slp_full.time, slptime_all{i});
    end
end
slp_full.pressure = double(slp_full.pressure);
slp_full.datetime = datetime(1800,01,01, 'TimeZone', 'UTC') + hours(slp_full.time); % converting to datetime
lati = double(ncread(nc_name, 'lat'));
loni = double(ncread(nc_name, 'lon'));
[slp_full.lat, slp_full.lon] = meshgrid(lati,loni); % transforms vectors lati and loni into arrays 
clear slpdata_all slptime_all i nc_name lati loni ii

end
