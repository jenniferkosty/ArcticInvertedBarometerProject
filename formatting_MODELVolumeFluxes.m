clc; clear;

%% Formatting model output data

%%% Strait Name
straits_all = {'BER', 'BAREN', 'FRAM', 'SPSIB'};

%%% Constructing array
start_day = datetime(2003, 1, 1, 'TimeZone', 'UTC');

%%% Set model run
u = 3;
label_options = {'Control', 'SLP', 'Wind'};
label = label_options{u};

%%% Loading hourly data, looping through the straits
for m = 1:length(straits_all)
    %%% Selecting strait
    strait_no = straits_all{m};

    %%% Loading data
    if strcmp(label, 'Control')
        file_label = 'CTR';
    elseif strcmp(label, 'SLP')
        file_label = 'SLP';
    elseif strcmp(label, 'Wind')
        file_label = 'WIND';
    end
    file_dir = 'C:\Users\jak279\OneDrive - Yale University\Research\Data\ModelOutput081126\';
    T = readtable(file_dir + string(strait_no) + '-2003-2023-ice-' + string(file_label) + '-Aug-6');

    %%% Removing lines with years embedded
    ind = T.Var1 > 2000;
    T = T(~ind,:);

    %%% Extracting volume flux and time data (time from the first strait)
    if m == 1
        volume_flux_full.datetime = NaT(numel(T), 1, 'TimeZone', 'UTC');
    end
    volume_flux_full.(string(strait_no)) = NaN(numel(T), 1);
    for i = 0:size(T, 1)-1

        start_idx = (i*24)+1;
        end_idx = start_idx + 23;
        if m == 1
            volume_flux_full.datetime(start_idx:end_idx) = start_day + days(i) + hours(0:23);
        end

        %%% Volume fluxes need to be multiplied by the cell width!
        cell_width = 55555;
        if strcmp(strait_no, 'BER')
            cell_width = 1*55555;
        end
        volume_flux_full.(string(strait_no))(start_idx:end_idx) = table2array(T(i+1,:)) .* cell_width;
    end
end

%%% Creating coarser time grid
dt = hours(6);
first_year = year(volume_flux_full.datetime(1));
last_year = year(volume_flux_full.datetime(end));
tgrid = datetime(first_year, 1, 1, "TimeZone", "UTC"):dt:datetime(last_year+1, 1, 1, "TimeZone", "UTC");
tgrid = tgrid(1:end-1);
tgrid = tgrid(:);

volume_flux.datetime = tgrid;

%%% Creating arrays to hold coarser volume flux data
for m = 1:length(straits_all)
    strait_no = straits_all{m};
    volume_flux.(string(strait_no)) = NaN(size(tgrid));
end

%%% Running mean, looping through all straits
for m = 1:length(straits_all)
    strait_no = straits_all{m};
    [~, volume_flux.(string(strait_no))] = running_mean(volume_flux_full.datetime, volume_flux_full.(string(strait_no)), tgrid);
end

dt_string = char(dt);
dt_string = string(dt_string(~isspace(dt_string)));
save('C:\Users\jak279\OneDrive - Yale University\Research\Data\OBP_interp\ModelOutput' + string(label) + 'VolumeFluxes_' + string(dt_string) + '.mat', 'volume_flux');



%%

figure()
hold on
for i = 1:length(straits_all)
    strait_no = straits_all{i};
    plot(volume_flux.datetime, movmean(volume_flux.(string(strait_no)), 4*30) ./ 1e6, 'DisplayName', string(upper(strait_no)));
end
legend()


