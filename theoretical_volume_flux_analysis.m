clc; clear;

%% Loading mooring data

moorings_all = {'a', 'b', 'd', 'np'};
dt_string = '1day';
dbar_to_pa = 10000;
for i = 1:length(moorings_all)

    %%% Selecting mooring
    mooring_no = moorings_all{i};

    %%% Loading BPR data
    data = load('C:\Users\jak279\OneDrive - Yale University\Research\Data\OBP_interp\BPR_mooring_' + string(mooring_no) + '_' + string(dt_string) + '.mat'); % BPR
    bpr.datetime = data.bpr.datetime;
    bpr_tmp = data.bpr.detided;
    bpr.(string(mooring_no) + '_detided') = bpr_tmp;
    bpr.(string(mooring_no) + '_start') = data.bpr.start;
    bpr.(string(mooring_no) + '_end') = data.bpr.end;
    bpr.(string(mooring_no) + '_lat') = data.bpr.lat;

    %%% Loading SLP data
    data = load('C:\Users\jak279\OneDrive - Yale University\Research\Data\OBP_interp\SLP_mooring_' + string(mooring_no) + '_' + string(dt_string) + '.mat'); % SLP
    slp.datetime = data.slp.datetime;
    slp.(string(mooring_no)) = NaN(size(data.slp.pressure));
    slp_tmp = data.slp.pressure .* dbar_to_pa; % kg/m/s^2
    for u = 1:length(bpr.(string(mooring_no) + '_start'))
        t1 = bpr.(string(mooring_no) + '_start')(u);
        t2 = bpr.(string(mooring_no) + '_end')(u);
        ind = data.slp.datetime >= t1 & data.slp.datetime <= t2;
        slp.(string(mooring_no))(ind) = slp_tmp(ind) - mean(slp_tmp(ind), 'omitnan');
    end
end

%% Loading bottom velocity data (and deriving fluxes)

load('C:\Users\jak279\OneDrive - Yale University\Research\Data\A3_daily_vel_1990to2024.mat')

A3_daily_vel.across = A3_daily_vel.u*cosd(315) - A3_daily_vel.v*sind(315); % across-strait
A3_daily_vel.across_vol = A3_daily_vel.across .* 4.25 .* 1000 .* 1000; % volume flux
A3_daily_vel.along = A3_daily_vel.u*sind(315) + A3_daily_vel.v*cosd(315); % along-strait
A3_daily_vel.along_vol = A3_daily_vel.along .* 4.25 .* 1000 .* 1000; % volume flux

A3_daily_vel.datetime = datetime(A3_daily_vel.time, 'ConvertFrom', 'datenum', 'TimeZone', 'UTC');

%% Computing required sea level height change (assuming perfect inverted barometer)

rho0 = 1027; %kg/m^3
g = 9.8; % m/s^2
for i = 1:length(moorings_all)
    mooring_no = moorings_all{i};
    slp.(string(mooring_no) + '_ib_ssh') = -(1./(rho0*g).*slp.(string(mooring_no)));
end

% %% Computing required volume flux assuming geostrophy
% 
% mooring_depth = [3825 3825 3530 4300];
% for i = 1:length(moorings_all)
%     mooring_no = moorings_all{i};
%     f = gsw_f(bpr.(string(mooring_no) + '_lat'));
%     slp.(string(mooring_no) + '_volume_flux') = -(g.*slp.(string(mooring_no) + '_ib_ssh').*mooring_depth(i)) ./ f;
% end

%% Computing required volume flux assuming mass conservation

mooring_no = 'a';
bg_area = 7e11;
arctic_area = 1.4e13;
time_scales = [3 10 30 100 365];
figure('Position', [10 10 700 700])
tiledlayout(3, 1)
lw = 2;

ax(1) = nexttile();
hold on
for i = 1:length(time_scales)
    volume_estimate_bg = (slp.(string(mooring_no) + '_ib_ssh') .* bg_area) ./ (time_scales(i)*24*60*60);
    plot(slp.datetime, volume_estimate_bg, 'LineWidth', lw, 'DisplayName', string(num2str(time_scales(i))) + ' days')

end
grid on
title('Canada Basin')
ylabel('Volume Flux [m^3/s]')
ylim(2e5.*[-1 1])
legend()

ax(2) = nexttile();
hold on
for i = 1:length(time_scales)
    volume_estimate_bg = (slp.(string(mooring_no) + '_ib_ssh') .* arctic_area) ./ (time_scales(i)*24*60*60);
    plot(slp.datetime, volume_estimate_bg, 'LineWidth', lw, 'DisplayName', string(num2str(time_scales(i))) + ' days')

end
grid on
title('Whole Arctic')
ylabel('Volume Flux [m^3/s]')
ylim(2e7.*[-1 1])

ax(3) = nexttile();
hold on
plot(A3_daily_vel.datetime, A3_daily_vel.along_vol - mean(A3_daily_vel.along_vol, 'omitnan'), 'LineWidth', lw, 'DisplayName', '"Observed" Volume Flux');
grid on
title('Observed')
ylabel('Volume Flux [m^3/s]')

linkaxes(ax, 'x')
xlim([datetime(2005, 1, 1, 'TimeZone', 'UTC') datetime(2006, 1, 1, 'TimeZone', 'UTC')])
ylim(4e8.*[-1 1])
%% Plotting data

figure('Position', [10 10 1000 500])
lw = 2;
yyaxis left
plot(slp.datetime, slp.a_volume_flux, 'LineWidth', lw, 'DisplayName', 'Required Volume Flux');
ylabel('Required Volume Flux for IB-Response [m^3/s]')
ylim(4e8.*[-1 1])
yyaxis right
plot(A3_daily_vel.datetime, A3_daily_vel.along_vol - mean(A3_daily_vel.along_vol, 'omitnan'), 'LineWidth', lw, 'DisplayName', '"Observed" Volume Flux');
xlim([datetime(2005, 1, 1, 'TimeZone', 'UTC') datetime(2006, 1, 1, 'TimeZone', 'UTC')])
ylabel('Volume Flux from Bering Strait Current Meter [m^3/s]')
grid on
ylim(4e8.*[-1 1])


%% Plotting spectra

figure()
spectral_estimate_pwelch(A3_daily_vel.datetime, A3_daily_vel.along_vol, 'Color','k', 'NumSegments', 8, 'MakeFigure', 1);
hold on
spectral_estimate_pwelch(slp.datetime, slp.a_volume_flux, 'Color', 'r', 'NumSegments', 8, 'MakeFigure', 1);
