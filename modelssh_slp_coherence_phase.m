clc; clear;

%% Loading model output data

moorings_all = {'a', 'b', 'd', 'np'};
dt_string = '6hr';
for i = 1:length(moorings_all)
    mooring_no = moorings_all{i};

    %%% Loading model output - Control
    data = load('C:\Users\jak279\OneDrive - Yale University\Research\Data\OBP_interp\ModelOutputControl_mooring_' + string(mooring_no) + '_' + string(dt_string) + '.mat');
    ssh_control.datetime = data.ssh.datetime;
    ssh_control.(string(mooring_no)) = data.ssh.ssh_anomaly ./ 100;
    ssh_control.(string(mooring_no) + '_start') = data.ssh.start;
    ssh_control.(string(mooring_no) + '_end') = data.ssh.end;

    %%% Loading model output - SLP only
    data = load('C:\Users\jak279\OneDrive - Yale University\Research\Data\OBP_interp\ModelOutputSLP_mooring_' + string(mooring_no) + '_' + string(dt_string) + '.mat');
    ssh_slp_only.datetime = data.ssh.datetime;
    ssh_slp_only.(string(mooring_no)) = data.ssh.ssh_anomaly ./ 100;
    ssh_slp_only.(string(mooring_no) + '_start') = data.ssh.start;
    ssh_slp_only.(string(mooring_no) + '_end') = data.ssh.end;

    %%% Loading model output - Wind only
    data = load('C:\Users\jak279\OneDrive - Yale University\Research\Data\OBP_interp\ModelOutputWind_mooring_' + string(mooring_no) + '_' + string(dt_string) + '.mat');
    ssh_wind_only.datetime = data.ssh.datetime;
    ssh_wind_only.(string(mooring_no)) = data.ssh.ssh_anomaly ./ 100;
    ssh_wind_only.(string(mooring_no) + '_start') = data.ssh.start;
    ssh_wind_only.(string(mooring_no) + '_end') = data.ssh.end;

    %%% Loading SLP
    data = load('C:\Users\jak279\OneDrive - Yale University\Research\Data\OBP_interp\ModelSLP_mooring_' + string(mooring_no) + '_' + string(dt_string) + '.mat');
    model_slp_forcing.datetime = data.slp.datetime;
    model_slp_forcing.(string(mooring_no)) = NaN(size(data.slp.pressure));
    pressure_tmp = data.slp.pressure;
    for u = 1:length(ssh_control.(string(mooring_no) + '_start'))
        t1 = ssh_control.(string(mooring_no) + '_start')(u);
        t2 = ssh_control.(string(mooring_no) + '_end')(u);
        ind = data.slp.datetime >= t1 & data.slp.datetime <= t2;
        model_slp_forcing.(string(mooring_no))(ind) = pressure_tmp(ind) - mean(pressure_tmp(ind), 'omitnan');
    end
    clear u t1 t2 ind pressure_tmp

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% VERSION 1: Computing coherence/transfer function/etc for each deployment and averaging %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Level of smoothing
P = 8;
numSegs = 8;
use_pwelch = 1;

no_runs = 3;

%%% Creating figure
fig = figure('Position', [10 10 900 800]);
tiledlayout(4, no_runs,'TileSpacing', 'compact', 'Padding','tight')
sgtitle('SSH_{MODEL} and SLP')
lw = 2; 
colors = orderedcolors("gem");

%%% Looping through 3 model runs
for uu = 1:3

if uu == 1
    ssh = ssh_control;
    label = 'Control';
elseif uu == 2
    ssh = ssh_slp_only;
    label = 'SLP Only';
elseif uu == 3
    ssh = ssh_wind_only;
    label = 'Wind Only';
end

%%% Creating frequency grid for averaging
a = 1/365;
if strcmp(dt_string, '6hr')
    b = 2;
elseif strcmp(dt_string, '1day')
    b = 1/2;
end
f = logspace(log10(a), log10(b), 100);

%%% Pre-allocating cells to hold the spectral analysis
[coherence_all, cospectrum_all, transfer_function_all, ...
    phase_all, coherence_crit_all,...
    cospectrum_crit_all, transfer_function_crit_all] = deal(cell(1, length(moorings_all)));
m = 1;
for i = 1:length(moorings_all)
    mooring_no = moorings_all{i};

    %%% Looping through each deployment
    for u = 1:length(ssh.(string(mooring_no) + '_start'))

        %%% Isolating annual segments
        t1 = ssh.(string(mooring_no) + '_start')(u);
        t2 = ssh.(string(mooring_no) + '_end')(u);
        idx = ssh.datetime >= t1 & ssh.datetime <= t2;
        idy = model_slp_forcing.datetime >= t1 & model_slp_forcing.datetime <= t2;
        t = ssh.datetime(idx);
        x = ssh.(string(mooring_no))(idx);
        y = model_slp_forcing.(string(mooring_no))(idy);
       
        %%% Skipping if segment is empty
        if isempty(x) | isempty(y)
            coherence_all{i}(u,:) = NaN(size(f));
            cospectrum_all{i}(u,:) = NaN(size(f));
            transfer_function_all{i}(u,:) = NaN(size(f));
            phase_all{i}(u,:) = NaN(size(f));
            coherence_crit_all{i}(u,:) = NaN(size(f));
            cospectrum_crit_all{i}(u,:) = NaN(size(f));
            transfer_function_crit_all{i}(u,:) = NaN(size(f));
            continue
        end

        %%% Computing coherence and phase
        if use_pwelch == 1
            [f_tmp, coherence_tmp, cospectrum_tmp, transfer_function_tmp, phase_tmp, ~, coherence_crit_tmp, cospectrum_crit_tmp, transfer_function_crit_tmp] = coherence_phase_pwelch(t, x, y,'NumSegments', numSegs, 'MonteCarloSignificance', 1);
            m = m + 1;
        else
            [f_tmp, coh_tmp, phase_tmp] = coherence_phase_jlab(t, x, y, 'P', P, 'PlotCoherence', 0);
        end

        %%% Skipping if segment was too short
        if isnan(f_tmp)
            coherence_all{i}(u,:) = NaN(size(f));
            cospectrum_all{i}(u,:) = NaN(size(f));
            transfer_function_all{i}(u,:) = NaN(size(f));
            phase_all{i}(u,:) = NaN(size(f));
            coherence_crit_all{i}(u,:) = NaN(size(f));
            cospectrum_crit_all{i}(u,:) = NaN(size(f));
            transfer_function_crit_all{i}(u,:) = NaN(size(f));
            continue
        end

        %%% Interpolating to common frequency grid for averaging
        idx = ~isnan(coherence_tmp);
        coherence_all{i}(u,:) = interp1(f_tmp(idx), coherence_tmp(idx), f);

        idx = ~isnan(cospectrum_tmp);
        cospectrum_all{i}(u,:) = interp1(f_tmp(idx), cospectrum_tmp(idx), f);

        idx = ~isnan(transfer_function_tmp);
        transfer_function_all{i}(u,:) = interp1(f_tmp(idx), transfer_function_tmp(idx), f);

        idx = ~isnan(phase_tmp);
        phase_all{i}(u,:) = interp1(f_tmp(idx), phase_tmp(idx), f);

        %%% Skipping if Monte Carlo was not conducted
        if isnan(coherence_crit_tmp)
            coherence_crit_all{i}(u,:) = NaN(size(f));
            cospectrum_crit_all{i}(u,:) = NaN(size(f));
            transfer_function_crit_all{i}(u,:) = NaN(size(f));
            continue
        end

        idx = ~isnan(coherence_crit_tmp);
        coherence_crit_all{i}(u,:) = interp1(f_tmp(idx), coherence_crit_tmp(idx), f);

        idx = ~isnan(cospectrum_crit_tmp);
        cospectrum_crit_all{i}(u,:) = interp1(f_tmp(idx), cospectrum_crit_tmp(idx), f);

        idx = ~isnan(transfer_function_crit_tmp);
        transfer_function_crit_all{i}(u,:) = interp1(f_tmp(idx), transfer_function_crit_tmp(idx), f);

    end
    
end


%%% Plotting coherence
ax = nexttile(uu);
title(label)
plot_mean_std_band(ax, f, coherence_all, moorings_all, colors, lw, ...
    'YLim',[0 1], 'YTicks',0:0.1:1, 'XLog',true, 'YLog',false, ...
    'YLabel', {'Magnitude-Squared', 'Coherence'}, ...
    'ShowLegend', 0, 'LegendLocation','northeast', ...
    'ShowXTicks', true, ...
    'TiePeriods',[100 20 3], 'TieLabels',{'100 days','20 days','3 days'}, 'ShowTieText', uu == 1);

%%% Plotting cross-spectrum
ax = nexttile(no_runs + uu);
plot_mean_std_band(ax, f, cospectrum_all, moorings_all, colors, lw, ...
    'YLim',[0 0.5], 'XLog',true, 'YLog',true, ...
    'YLabel', "Cross-Spectrum", ...
    'ShowXTicks', false, 'TiePeriods',[100 20 3], 'ShowTieText', false);

%%% Plotting transfer function
ax = nexttile(2*no_runs + uu);
plot_mean_std_band(ax, f, transfer_function_all, moorings_all, colors, lw, ...
    'YLim',[-1 1], 'YTicks',-1:0.25:1, 'XLog',true, 'YLog',false, ...
    'YLabel',  "Transfer Function", ...
    'ShowXTicks', false, 'TiePeriods',[100 20 3], 'ShowTieText', false);
yline(ax, 0, '-', [0.5 0.5 0.5], 'LineWidth', lw);

%%% Plotting phase
ax = nexttile(3*no_runs+uu);
phase_all_corrected = phase_all;
for i = 1:length(moorings_all)
    ind = coherence_all{i} < coherence_crit_all{i};
    phase_all_corrected{i}(ind) = NaN;
end
plot_mean_std_band(ax, f, phase_all_corrected, moorings_all, colors, lw, ...
    'YLim',[-180 180], 'YTicks',-180:60:180, 'XLog',true, 'YLog',false, ...
    'XLabel',"Frequency [cycles/day]", ...
    'YLabel', "Phase [deg]", ...
    'ShowXTicks', true, 'TiePeriods',[100 20 3], 'ShowTieText', false, 'isPhase', 1, 'PhaseBandPrct', [16 84]);
yline(ax, 0, '-', [0.5 0.5 0.5], 'LineWidth', lw);

end

%%% Saving figure
if use_pwelch == 1
    print(fig,'C:\Users\jak279\OneDrive - Yale University\Research\Figures\inverted_barometer_paper\modeloutput_slp_coherence_phase_pwelch.png','-dpng', '-r600')
else 
    print(fig,'C:\Users\jak279\OneDrive - Yale University\Research\Figures\inverted_barometer_paper\modeloutput_slp_coherence_phase_jlab.png','-dpng', '-r600')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Version 2: Computing average spectra for all years and then computing coherence %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numSegs = 8;
numSurr = 50;

no_runs = 3;

%%% Creating figure
fig = figure('Position', [10 10 900 800]);
tiledlayout(3, no_runs,'TileSpacing', 'compact', 'Padding','tight')
sgtitle('SSH_{MODEL} and SLP')
lw = 2; fs = 12;
colors = orderedcolors("gem");

%%% Looping through 3 model runs
for uu = 1:no_runs

if uu == 1
    ssh = ssh_control;
    label = 'Control';
elseif uu == 2
    ssh = ssh_slp_only;
    label = 'SLP Only';
elseif uu == 3
    ssh = ssh_wind_only;
    label = 'Wind Only';
end

%%% Pre-allocating cells to hold the spectral analysis
[coherence_all, cospectrum_all, transfer_function_all, ...
    phase_all, theor_coherence_crit_all, coherence_crit_all,...
    cospectrum_crit_all, transfer_function_crit_all] = deal(cell(1, length(moorings_all)));

%%% Looping through moorings
for i = 1:length(moorings_all)
    mooring_no = moorings_all{i};
    ind = model_slp_forcing.datetime >= ssh.datetime(1) & model_slp_forcing.datetime <= ssh.datetime(end);
    tx = model_slp_forcing.datetime(ind);
    x = model_slp_forcing.(string(mooring_no))(ind);
    ty = ssh.datetime;
    y = ssh.(string(mooring_no));
    [f, coherence_all{i}, cospectrum_all{i}, transfer_function_all{i}, phase_all{i}, theor_coherence_crit_all{i}, coherence_crit_all{i}, cospectrum_crit_all{i}, transfer_function_crit_all{i}] = coherence_phase_pwelch_full_record(tx, x, ty, y, 'StartDates', ssh.(string(mooring_no) + '_start'), 'EndDates', ssh.(string(mooring_no) + '_end'), 'NumSegments', numSegs, 'Nsurr', numSurr, 'MonteCarloSignificance', 1);
end

%%% Plotting coherence
ax = nexttile(0*no_runs+uu);
plot_mean_sig_thresh(ax, f, coherence_all, coherence_crit_all, moorings_all, colors, lw, ...
    'YLim',[0 1], 'YTicks',0:0.25:1, 'XLog',true, 'YLog',false, ...
    'YLabel', {'Magnitude-Squared', 'Coherence'}, ...
    'ShowLegend', 0, 'LegendLocation','northeast', ...
    'ShowXTicks', false, ...
    'TiePeriods',[100 20 3], 'TieLabels',{'100 days','20 days','3 days'}, 'ShowTieText', true, 'Title', label);
if uu ~= 1
    set(ax, 'YTickLabel', []);
    ylabel(ax, '');
end

%%% Transfer Function
ax = nexttile(1*no_runs+uu);
plot_mean_sig_thresh(ax, f, transfer_function_all, transfer_function_crit_all, moorings_all, colors, lw, ...
    'YLim',1.15.*[0 1], 'YTicks',-1:0.25:1, 'XLog',true, 'YLog',false, ...
    'YLabel',  "Transfer Function", ...
    'ShowXTicks', false, 'TiePeriods',[100 20 3], 'ShowTieText', false);
yline(ax, 0, '-', 'Color', [0.5 0.5 0.5], 'LineWidth', lw);
if uu ~= 1
    set(ax, 'YTickLabel', []);
    ylabel(ax, '');
end

%%% Phase
ax = nexttile(2*no_runs+uu);
phase_all_corrected = phase_all;
for i = 1:length(moorings_all)
    ind = coherence_all{i} < coherence_crit_all{i};
    phase_all_corrected{i}(ind) = NaN;
    phase_all_corrected{i} = rad2deg(unwrap(phase_all_corrected{i}));
end
plot_mean_sig_thresh(ax, f, phase_all_corrected, cell(size(phase_all_corrected)), moorings_all, colors, lw, ...
    'YLim',[-180 180], 'YTicks',-180:60:180, 'XLog',true, 'YLog',false, ...
    'XLabel',"Frequency [cycles/day]", ...
    'YLabel', "Phase [deg]", ...
    'ShowXTicks', true, 'TiePeriods',[100 20 3], 'ShowTieText', false);
yline(ax, 0, '-', 'Color', [0.5 0.5 0.5], 'LineWidth', lw);
if uu ~= 1
    set(ax, 'YTickLabel', []);
    ylabel(ax, '');
end

end

print(fig,'C:\Users\jak279\OneDrive - Yale University\Research\Figures\inverted_barometer_paper\modeloutput_slp_coherence_phase_pwelch.png','-dpng', '-r600')

%% Bandpassing data

minutes_per_hour = 60;
hours_per_day = 24;
cutoff_lower = hours_per_day * 3; % in hours
cutoff_upper = hours_per_day * 20; % in hours

%%% SSH Control at mooring sites
dt = hours(ssh_control.datetime(2) - ssh_control.datetime(1));
sampling_interval = minutes_per_hour * dt; % in minutes
for i = 1:length(moorings_all)
    mooring_no = moorings_all(i);
    [ssh_control.(string(mooring_no) + '_highpassed'), ~] = highpassing_butterworth(ssh_control.(string(mooring_no)), sampling_interval, cutoff_lower);
    [tmp, ssh_control.(string(mooring_no) + '_lowpassed')] = highpassing_butterworth(ssh_control.(string(mooring_no)), sampling_interval, cutoff_upper);
    [~, ssh_control.(string(mooring_no) + '_bandpassed')] = highpassing_butterworth(tmp, sampling_interval, cutoff_lower);
end

%%% SSH SLP Only at mooring sites
dt = hours(ssh_slp_only.datetime(2) - ssh_slp_only.datetime(1));
sampling_interval = minutes_per_hour * dt; % in minutes
for i = 1:length(moorings_all)
    mooring_no = moorings_all(i);
    [ssh_slp_only.(string(mooring_no) + '_highpassed'), ~] = highpassing_butterworth(ssh_slp_only.(string(mooring_no)), sampling_interval, cutoff_lower);
    [tmp, ssh_slp_only.(string(mooring_no) + '_lowpassed')] = highpassing_butterworth(ssh_slp_only.(string(mooring_no)), sampling_interval, cutoff_upper);
    [~, ssh_slp_only.(string(mooring_no) + '_bandpassed')] = highpassing_butterworth(tmp, sampling_interval, cutoff_lower);
end


%%% SSH Wind Only at mooring sites
dt = hours(ssh_wind_only.datetime(2) - ssh_wind_only.datetime(1));
sampling_interval = minutes_per_hour * dt; % in minutes
for i = 1:length(moorings_all)
    mooring_no = moorings_all(i);
    [ssh_wind_only.(string(mooring_no) + '_highpassed'), ~] = highpassing_butterworth(ssh_wind_only.(string(mooring_no)), sampling_interval, cutoff_lower);
    [tmp, ssh_wind_only.(string(mooring_no) + '_lowpassed')] = highpassing_butterworth(ssh_wind_only.(string(mooring_no)), sampling_interval, cutoff_upper);
    [~, ssh_wind_only.(string(mooring_no) + '_bandpassed')] = highpassing_butterworth(tmp, sampling_interval, cutoff_lower);
end


%%% SLP at mooring sites
dt = hours(model_slp_forcing.datetime(2) - model_slp_forcing.datetime(1));
sampling_interval = minutes_per_hour * dt; % in minutes
for i = 1:length(moorings_all)
    mooring_no = moorings_all(i);
    [model_slp_forcing.(string(mooring_no) + '_highpassed'), ~] = highpassing_butterworth(model_slp_forcing.(string(mooring_no)), sampling_interval, cutoff_lower);
    [tmp, model_slp_forcing.(string(mooring_no) + '_lowpassed')] = highpassing_butterworth(model_slp_forcing.(string(mooring_no)), sampling_interval, cutoff_upper);
    [~, model_slp_forcing.(string(mooring_no) + '_bandpassed')] = highpassing_butterworth(tmp, sampling_interval, cutoff_lower);
end

%%

figure()
tiledlayout(2, 1)

mooring_no = 'a';
tx = model_slp_forcing.datetime;
x = model_slp_forcing.(string(mooring_no));
ty = ssh_control.datetime;
y = ssh_control.(string(mooring_no));

%%% Wavelet transfer
ax(1) = nexttile(1);
ax(2) = nexttile(2);
plot_wavelet_transfer(ax(1), ax(2), tx, x, y, 'PlotPhase', true, 'RunMonteCarlo', false, ...
    'NSurrogates', 50, 'SigLevel', 95, 'YTicks', [1 4 16 64 256 1024], 'LogGain', false);
hold on
clim(ax(1), [0 1])
%xlim(ax, [2003 2010])

%%
fig = figure('Position', [10 10 1200 600]);
tiledlayout(2, 1, 'Padding', 'tight', 'TileSpacing', 'tight')
clear ax
mooring_no = 'a';
ind = 11;
t1 = ssh_control.(string(mooring_no) + '_start')(ind);
t2 = ssh_control.(string(mooring_no) + '_end')(ind);
c1 = [1 0.6 0];


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Wavelet Cross-Power %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

tx = ssh_control.datetime;
x = ssh_control.(string(mooring_no));
ty = model_slp_forcing.datetime;
y = model_slp_forcing.(string(mooring_no));

%%% Wavelet coherence
ax(1) = nexttile();
plot_wavelet_coherence(ax(1), [], tx, x, y, 'PlotPhase', false, 'RunMonteCarlo', true, ...
    'NSurrogates', 50, 'SigLevel', 95, 'YTicks', [1 4 16 64 256 1024]);
hold on

xline(yearfrac(datenum([t1 t2])), 'Color', c1, 'LineWidth', lw)
xlabel(' ')
set(gca, 'FontSize', fs)
text(0.02, 0.93, '(a)', 'Units', 'normalized', 'Color', 'r',...
    'HorizontalAlignment', 'center', 'FontSize', fs, ...
    'FontWeight', 'bold', 'BackgroundColor', 'none', 'Parent', ax(1))


%%%%%%%%%%%%%%%%%%%%%%%
%%% Bandpassed Data %%%
%%%%%%%%%%%%%%%%%%%%%%%

ax(2) = nexttile();
yyaxis left
hold on
plot(ssh_control.datetime, -ssh_control.(string(mooring_no) + '_bandpassed'), 'LineWidth', lw)
ylabel({'-SSH Anomaly', '[dbar]'})
ylim(0.3.*[-1 1])
yyaxis right
hold on
plot(model_slp_forcing.datetime, model_slp_forcing.(string(mooring_no) + '_bandpassed'), 'LineWidth', lw)
ylabel({'SLP Anomaly', '[dbar]'}, 'Rotation', 270)
ylim(0.3.*[-1 1])
grid on
xlim([t1 t2])
set(gca, 'FontSize', fs)
yline(0, 'Color', [0.5 0.5 0.5], 'LineWidth', lw)
text(0.5, 0.93, 'Bandpassed 3-20 days', 'Units', 'normalized', ...
    'HorizontalAlignment', 'center', 'FontSize', fs, ...
    'FontWeight', 'bold', 'BackgroundColor', 'none', 'Parent', ax(2))
text(0.02, 0.93, '(b)', 'Units', 'normalized', 'Color', 'r',...
    'HorizontalAlignment', 'center', 'FontSize', fs, ...
    'FontWeight', 'bold', 'BackgroundColor', 'none', 'Parent', ax(2))

% Get positions of the axes in figure coordinates [x y width height]
pos_slp  = ax(1).Position;   % BPR-SLP cross-wavelet
pos_ts   = ax(2).Position;   % Time series

% Get the x-limits of the BPR-SLP cross-wavelet to convert t1,t2 to normalized coords
xlim_slp = ax(1).XLim;
t1_norm  = yearfrac(datenum(t1));
t2_norm  = yearfrac(datenum(t2));

% Convert t1 and t2 to figure x coordinates within the SLP wavelet axes
x1_fig = pos_slp(1) + (t1_norm - xlim_slp(1)) / (xlim_slp(2) - xlim_slp(1)) * pos_slp(3);
x2_fig = pos_slp(1) + (t2_norm - xlim_slp(1)) / (xlim_slp(2) - xlim_slp(1)) * pos_slp(3);

% Bottom of SLP wavelet axes
y_slp_bottom = pos_slp(2);

% Top of time series axes
y_ts_top = pos_ts(2) + pos_ts(4);

% Left edge of time series axes (t1)
x1_ts_fig = pos_ts(1);
% Right edge of time series axes (t2)
x2_ts_fig = pos_ts(1) + pos_ts(3);

% Create invisible overlay axes covering the whole figure
ax_overlay = axes('Position', [0 0 1 1], 'Visible', 'off');
ax_overlay.XLim = [0 1];
ax_overlay.YLim = [0 1];
hold(ax_overlay, 'on')

% Draw the connecting lines with alpha
l1 = plot(ax_overlay, [x1_fig, x1_ts_fig], [y_slp_bottom, y_ts_top], ...
    '-', 'Color', c1, 'LineWidth', lw);
l2 = plot(ax_overlay, [x2_fig, x2_ts_fig], [y_slp_bottom, y_ts_top], ...
    '-', 'Color', c1, 'LineWidth', lw);

% Set alpha
l1.Color(4) = 0.4;   % 0 = fully transparent, 1 = fully opaque
l2.Color(4) = 0.4;
