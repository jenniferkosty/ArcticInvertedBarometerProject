%% Script to format Arctic volume transport estimates

clc; clear;

%% Loading data

directory = 'C:\Users\jak279\OneDrive - Yale University\Research\Data\Tsubuoschi_all\';
filename = string(directory) + 'ocean_VFacs_Oct2004toMay2010_170123_VabsScons_var_mon.txt';

%%% Open and detect where data starts (skip comment/header lines)
fid = fopen(filename, 'r');
headerLines = 0;
while true
    line = fgetl(fid);
    if ischar(line) && ~isempty(line) && ~isempty(regexp(line, '^\s*\d{4}', 'once'))
        break;  % Found the first data line (starts with a 4-digit year)
    end
    headerLines = headerLines + 1;
end
fclose(fid);

%%% Importing data after the header
data = readmatrix(filename, 'NumHeaderLines', headerLines, 'Delimiter', ',');

%%% Formatting in a structure
strait.year    = data(:, 1);
strait.month   = data(:, 2);
strait.net     = data(:, 3);
strait.davis   = data(:, 4);
strait.fram    = data(:, 5);
strait.bso     = data(:, 6);
strait.bering  = data(:, 7);
strait.belgica = data(:, 8);
strait.egc     = data(:, 9);
strait.middle  = data(:, 10);
strait.wsc     = data(:, 11);

%%% Adding datetime array
strait.datetime = datetime(strait.year, strait.month, 1, 'TimeZone', 'UTC');

%%% Adding metadata
strait.units    = 'Sv  (1 Sv = 10^6 m^3 s^-1)';
strait.note     = 'Positive = inflow, Negative = outflow';
strait.source   = filename;

%%% Saving data
save('C:\Users\jak279\OneDrive - Yale University\Research\Data\OBP_interp\VolumeTransportEstimates.mat', 'strait');


%%

figure()
plot(strait.datetime, strait.bering)