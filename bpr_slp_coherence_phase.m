%% Script to visualize coherence/phase between OBP and SLP

clc; clear;

%% Loading mooring data

moorings_all = {'a', 'b', 'd', 'np'};
dt_string = '6hr';
use_ssh = 0;

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

    %%% Loading SLP data
    data = load('C:\Users\jak279\OneDrive - Yale University\Research\Data\OBP_interp\SLP_mooring_' + string(mooring_no) + '_' + string(dt_string) + '.mat'); % SLP
    slp.datetime = data.slp.datetime;
    slp_tmp = data.slp.pressure;
    %slp.(string(mooring_no)) = slp_tmp;
    %slp_tmp = -gsw_z_from_p(slp_tmp, data.slp.actual_lat);
    slp.(string(mooring_no)) = NaN(size(slp_tmp));
    for u = 1:length(bpr.(string(mooring_no) + '_start'))
        t1 = bpr.(string(mooring_no) + '_start')(u);
        t2 = bpr.(string(mooring_no) + '_end')(u);
        ind = data.slp.datetime >= t1 & data.slp.datetime <= t2;
        slp.(string(mooring_no))(ind) = slp_tmp(ind) - mean(slp_tmp(ind), 'omitnan');
    end
end

%%% Computing SSH_b = OBP - SLP
for i = 1:length(moorings_all)
    mooring_no = moorings_all{i};
    bpr.(string(mooring_no) + '_ssh') = bpr.(string(mooring_no) + '_detided') - slp.(string(mooring_no));
end

clear bpr_tmp data i ind slp_tmp t1 t2 u 

%%
figure()
tiledlayout(2, 1)

mooring_no = 'np';
tx = slp.datetime;
x = slp.(string(mooring_no));
ty = bpr.datetime;
y = bpr.(string(mooring_no) + '_detided');

%%% Wavelet transfer
ax(1) = nexttile(1);
ax(2) = nexttile(2);
plot_wavelet_transfer(ax(1), ax(2), ty, x, y, 'PlotPhase', true, 'RunMonteCarlo', false, ...
    'NSurrogates', 50, 'SigLevel', 95, 'YTicks', [1 4 16 64 256 1024], 'LogGain', false);
hold on
clim(ax(1), [0 1])
%xlim(ax, [2003 2010])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% VERSION 1: Computing coherence/transfer function/etc for each deployment and averaging %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Level of smoothing
P = 8;
numSegs = 8;
use_pwelch = 1;

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

for i = 1:length(moorings_all)
    mooring_no = moorings_all{i};
    
    %%% Looping through each deployent period
    for u = 1:length(bpr.(string(mooring_no) + '_start'))

        %%% Isolating annual segments
        t1 = bpr.(string(mooring_no) + '_start')(u);
        t2 = bpr.(string(mooring_no) + '_end')(u);
        idx = bpr.datetime >= t1 & bpr.datetime <= t2;
        t = bpr.datetime(idx);
        x = bpr.(string(mooring_no) + '_detided')(idx);
        var1 = 'obp';
        idy = slp.datetime >= t1 & slp.datetime <= t2;
        y = slp.(string(mooring_no))(idy);
        var2 = 'slp';
       

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
        else
            [f_tmp, coherence_tmp, cospectrum_tmp, transfer_function_tmp, phase_tmp, coherence_crit_tmp, cospectrum_crit_tmp, transfer_function_crit_tmp] = coherence_phase_jlab(t, x, y,'P', P, 'MonteCarloSignificance', 1);
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

%%% Plotting
figure('Position', [10 10 1000 500])
tiledlayout(1, 3, 'TileSpacing', 'tight', 'Padding', 'none');
colors = orderedcolors('gem');
lw = 2;

%%% Coherence
ax = nexttile(1);
plot_mean_std_band(ax, f, coherence_all, moorings_all, colors, lw, ...
    'YLim',[0 1], 'YTicks',0:0.1:1, 'XLog',true, 'YLog',false, ...
    'YLabel', 'Magnitude-Squared Coherence', ...
    'ShowLegend', 0, 'LegendLocation','northeast', ...
    'ShowXTicks', true, ...
    'TiePeriods',[100 20 3], 'TieLabels',{'100 days','20 days','3 days'}, 'ShowTieText', true);

%%% Transfer Function
ax = nexttile(2);
plot_mean_std_band(ax, f, transfer_function_all, moorings_all, colors, lw, ...
    'YLim',[-1 1], 'YTicks',-1:0.25:1, 'XLog',true, 'YLog',false, ...
    'YLabel',  "Transfer Function", ...
    'ShowXTicks', false, 'TiePeriods',[100 20 3], 'ShowTieText', false);
yline(ax, 0, '-', [0.5 0.5 0.5], 'LineWidth', lw);

%%% Phase
ax = nexttile(3);
phase_all_corrected = phase_all;
for i = 1:length(moorings_all)
    ind = coherence_all{i} < coherence_crit_all{i};
    phase_all_corrected{i}(ind) = NaN;
end
plot_mean_std_band(ax, f, phase_all_corrected, moorings_all, colors, lw, ...
    'YLim',[-180 180], 'YTicks',-180:60:180, 'XLog',true, 'YLog',false, ...
    'XLabel',"Frequency [cycles/day]", ...
    'YLabel', "Phase [deg]", ...
    'ShowXTicks', true, 'TiePeriods',[100 20 3], 'ShowTieText', false, 'isPhase', 1);
yline(ax, 0, '-', [0.5 0.5 0.5], 'LineWidth', lw);

print(fig,'C:\Users\jak279\OneDrive - Yale University\Research\Figures\inverted_barometer_paper\' + string(var1) + '_' + string(var2) + '_coherence_phase_pwelch.png','-dpng', '-r600')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Version 2: Computing average spectra for all years and then computing coherence %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numSegs = 8;
numSurr = 50;

%%% Pre-allocating cells to hold the spectral analysis
[coherence_all, cospectrum_all, transfer_function_all, ...
    phase_all, theor_coherence_crit_all, coherence_crit_all,...
    cospectrum_crit_all, transfer_function_crit_all] = deal(cell(1, length(moorings_all)));

%%% Looping through moorings
for i = 1:length(moorings_all)
    mooring_no = moorings_all{i};
    tx = slp.datetime; 
    x = slp.(string(mooring_no));
    var1 = 'slp';
    ty = bpr.datetime;
    y = bpr.(string(mooring_no) + '_detided');
    var2 = 'obp';
    [f, coherence_all{i}, cospectrum_all{i}, transfer_function_all{i}, phase_all{i}, theor_coherence_crit_all{i}, coherence_crit_all{i}, cospectrum_crit_all{i}, transfer_function_crit_all{i}] = coherence_phase_pwelch_full_record(tx, x, ty, y, 'StartDates', bpr.(string(mooring_no) + '_start'), 'EndDates', bpr.(string(mooring_no) + '_end'), 'NumSegments', numSegs, 'Nsurr', numSurr, 'MonteCarloSignificance', 1);
end
clear tx x ty y

%%% Plotting
fig = figure('Position', [10 10 1000 400]);
tiledlayout(1, 3, 'TileSpacing', 'tight', 'Padding', 'none');
%sgtitle('Relationship between OBP at mooring sites and OBP at Mooring D')
colors = orderedcolors('gem');
lw = 2;

%%% Coherence
ax = nexttile();
plot_mean_sig_thresh(ax, f, coherence_all, coherence_crit_all, moorings_all, colors, lw, ...
    'YLim',[0 1], 'YTicks',0:0.2:1, 'XLog',true, 'YLog',false, ...
    'YLabel', 'Magnitude-Squared Coherence', ...
    'XLabel', 'Frequency [cycles/day]', ...
    'Title', '(a)', ...
    'ShowLegend', 0, 'LegendLocation','southwest', ...
    'ShowXTicks', true, ...
    'TiePeriods',[100 20 3], 'TieLabels',{'100 days','20 days','3 days'}, 'ShowTieText', true);

%%% Transfer Function
ax = nexttile();
plot_mean_sig_thresh(ax, f, transfer_function_all, transfer_function_crit_all, moorings_all, colors, lw, ...
    'YLim',[0 1.1], 'YTicks',-1:0.2:1, 'XLog',true, 'YLog',false, ...
    'YLabel',  "Transfer Function Magnitude", ...
    'XLabel', "Frequency [cycles/day]", ...
    'Title', '(b)', ...
    'ShowXTicks', true, 'TiePeriods',[100 20 3], 'ShowTieText', false);

%%% Phase
ax = nexttile();
phase_all_corrected = phase_all;
for i = 1:size(phase_all_corrected, 2)
    ind = coherence_all{i} < coherence_crit_all{i};
    phase_all_corrected{i}(ind) = NaN;
    phase_all_corrected{i} = rad2deg(unwrap(phase_all_corrected{i}));
end
plot_mean_sig_thresh(ax, f, phase_all_corrected, cell(size(phase_all_corrected)), moorings_all, colors, lw, ...
    'YLim',[-180 180], 'YTicks',-180:60:180, 'XLog',true, 'YLog',false, ...
    'XLabel',"Frequency [cycles/day]", ...
    'YLabel', "Phase [deg]", ...
    'Title', '(c)', ...
    'ShowLegend', 1, 'LegendLocation','southeast', ...
    'ShowXTicks', true, 'TiePeriods',[100 20 3], 'ShowTieText', false);
yline(ax, 0, '-', 'Color', [0.5 0.5 0.5], 'LineWidth', lw, 'HandleVisibility', 'off');

%%% Saving figure
print(fig,'C:\Users\jak279\OneDrive - Yale University\Research\Figures\inverted_barometer_paper\' + string(var1) + '_' + string(var2) + '_coherence_phase_pwelch.png','-dpng', '-r600')

%% Plotting wavelet transforms

figure('Position', [10 10 800 700])
tiledlayout(2, 1)

ax(1) = nexttile();
x = bpr.a_detided;
t = bpr.datetime;
dt = datenum(t(2)) - datenum(t(1));
fs_num = 1./dt;
plot_wavelet(ax(1), t, x, fs_num, 'Title', 'OBP', 'CLim', log10([0 4]));

ax(2) = nexttile();
x = slp.a;
t = slp.datetime;
dt = datenum(t(2)) - datenum(t(1));
fs_num = 1./dt;
plot_wavelet(ax(2), t, x, fs_num, 'Title', 'SLP', 'Clim', log10([0 4]));

linkaxes(ax, 'x')
xlim([2003.5 2013.5])

%% Plotting wavelet coherence

fig = figure('Position', [10 10 800 700]);
tiledlayout(2,1, 'TileSpacing', 'tight', 'Padding', 'none')
mooring_no = 'a';
sgtitle('Relationship between OBP (x) and SLP (y) at Mooring ' + string(upper(mooring_no)))
var1 = 'bprA';
var2 = 'slpA';

%%% Isolating data
ind = slp.datetime >= bpr.datetime(1) & slp.datetime <= bpr.datetime(end);
tx = bpr.datetime;
x = bpr.(string(mooring_no) + '_detided');
ty = slp.datetime(ind);
y = slp.(string(mooring_no))(ind);

ax1 = nexttile(1);
ax2 = nexttile(2);
plot_wavelet_coherence(ax1, ax2, tx, x, y, 'RunMonteCarlo', true, 'NSurrogates', 50, 'SigLevel',95);
xlabel(ax1, ' ')
set(ax1, 'XTickLabel', ' ');
linkaxes([ax1 ax2], 'x')
%xlim([2003.5 2013.5])

%%% Saving figure
print(fig,'C:\Users\jak279\OneDrive - Yale University\Research\Figures\inverted_barometer_paper\' + string(var1) + '_' + string(var2) + '_wavelet_coherence.png','-dpng', '-r600')

%% Bandpassing data

minutes_per_hour = 60;
hours_per_day = 24;
cutoff_lower = hours_per_day * 3; % in hours
cutoff_upper = hours_per_day * 20; % in hours

%%% OBP at mooring sites
dt = hours(bpr.datetime(2) - bpr.datetime(1));
sampling_interval = minutes_per_hour * dt; % in minutes
for i = 1:length(moorings_all)
    mooring_no = moorings_all(i);
    [bpr.(string(mooring_no) + '_highpassed'), ~] = highpassing_butterworth(bpr.(string(mooring_no) + '_detided'), sampling_interval, cutoff_lower);
    [tmp, bpr.(string(mooring_no) + '_lowpassed')] = highpassing_butterworth(bpr.(string(mooring_no) + '_detided'), sampling_interval, cutoff_upper);
    [~, bpr.(string(mooring_no) + '_bandpassed')] = highpassing_butterworth(tmp, sampling_interval, cutoff_lower);
end

%%% SLP at mooring sites
dt = hours(slp.datetime(2) - slp.datetime(1));
sampling_interval = minutes_per_hour * dt; % in minutes
for i = 1:length(moorings_all)
    mooring_no = moorings_all(i);
    [slp.(string(mooring_no) + '_highpassed'), ~] = highpassing_butterworth(slp.(string(mooring_no)), sampling_interval, cutoff_lower);
    [tmp, slp.(string(mooring_no) + '_lowpassed')] = highpassing_butterworth(slp.(string(mooring_no)), sampling_interval, cutoff_upper);
    [~, slp.(string(mooring_no) + '_bandpassed')] = highpassing_butterworth(tmp, sampling_interval, cutoff_lower);
end

%% Time series showing relationship between OBP/SLP

mooring_no = 'a';
ind = 2;
t1 = bpr.(string(mooring_no) + '_start')(ind);
t2 = bpr.(string(mooring_no) + '_end')(ind);

figure('Position', [10 10 1000 700])

% Define positions [left bottom width height]
ax1 = subplot('Position', [0.1  0.73  0.8  0.21]);
ax2 = subplot('Position', [0.1  0.52  0.8  0.21]);
ax3 = subplot('Position', [0.1  0.31  0.8  0.21]);
ax4 = subplot('Position', [0.1  0.05  0.85  0.19]);

lw = 2;
fs = 12;

%%% Plotting time series
axes(ax1)
yyaxis left
hold on
plot(bpr.datetime, bpr.(string(mooring_no) + '_detided'), 'LineWidth', lw, 'DisplayName', 'OBP')
ylabel({'OBP Anomaly', '[dbar]'})
ylim(0.2.*[-1 1])
yticks(-0.1:0.1:0.1)
yyaxis right
hold on
plot(slp.datetime, slp.(string(mooring_no)), 'LineWidth', lw, 'DisplayName', 'SLP')
ylabel({'SLP Anomaly', '[dbar]'}, 'Rotation', 270)
ylim(0.4.*[-1 1])
yticks(-0.2:0.2:0.2)
grid on
xlim([t1 t2])
set(gca, 'FontSize', fs)
yline(0, 'Color', [0.5 0.5 0.5], 'LineWidth', lw)
xticklabels(' ')
text(0.02, 0.93, '(a)', 'Units', 'normalized', 'Color', 'r',...
    'HorizontalAlignment', 'center', 'FontSize', fs, ...
    'FontWeight', 'bold', 'BackgroundColor', 'none', 'Parent', ax1)

%%% Seasonal Cycle
axes(ax2)
yyaxis left
hold on
ind = bpr.datetime >= t1 & bpr.datetime <= t2;
x = bpr.(string(mooring_no) + '_detided')(ind);
t = bpr.datetime(ind);
cyc = [1 2 3];
harmonics = computing_harmonics(cyc, x, t);
plot(harmonics.time, harmonics.reconstruction, 'LineWidth', lw)
ylabel({'OBP Anomaly', '[dbar]'})
ylim(0.1.*[-1 1])
yyaxis right
hold on
x = slp.(string(mooring_no))(ind);
t = slp.datetime(ind);
cyc = [1 2 3];
harmonics = computing_harmonics(cyc, x, t);
plot(harmonics.time, harmonics.reconstruction, 'LineWidth', lw)
ylabel({'SLP Anomaly', '[dbar]'}, 'Rotation', 270)
ylim(0.1.*[-1 1])
grid on
xlim([t1 t2])
set(gca, 'FontSize', fs)
yline(0, 'Color', [0.5 0.5 0.5], 'LineWidth', lw)
xticklabels(' ')
text(0.5, 0.93, 'Seasonal Cycle', 'Units', 'normalized', ...
    'HorizontalAlignment', 'center', 'FontSize', fs, ...
    'FontWeight', 'bold', 'BackgroundColor', 'none', 'Parent', ax2)
text(0.02, 0.93, '(b)', 'Units', 'normalized', 'Color', 'r',...
    'HorizontalAlignment', 'center', 'FontSize', fs, ...
    'FontWeight', 'bold', 'BackgroundColor', 'none', 'Parent', ax2)

%%% Bandpassed
axes(ax3);
yyaxis left
hold on
plot(bpr.datetime, bpr.(string(mooring_no) + '_bandpassed'), 'LineWidth', lw)
ylabel({'OBP Anomaly', '[dbar]'})
ylim(0.15.*[-1 1])
yyaxis right
hold on
plot(slp.datetime, slp.(string(mooring_no) + '_bandpassed'), 'LineWidth', lw)
ylabel({'SLP Anomaly', '[dbar]'}, 'Rotation', 270)
ylim(0.3.*[-1 1])
grid on
xlim([t1 t2])
set(gca, 'FontSize', fs)
yline(0, 'Color', [0.5 0.5 0.5], 'LineWidth', lw)
% *** No xticklabels(' ') here so x-axis labels are visible ***
text(0.5, 0.93, 'Bandpassed 3-20 days', 'Units', 'normalized', ...
    'HorizontalAlignment', 'center', 'FontSize', fs, ...
    'FontWeight', 'bold', 'BackgroundColor', 'none', 'Parent', ax3)
text(0.02, 0.93, '(c)', 'Units', 'normalized', 'Color', 'r',...
    'HorizontalAlignment', 'center', 'FontSize', fs, ...
    'FontWeight', 'bold', 'BackgroundColor', 'none', 'Parent', ax3)

%%% Formatting data for wavelet coherence
ind = slp.datetime >= bpr.datetime(1) & slp.datetime <= bpr.datetime(end);
tx = bpr.datetime;
x = bpr.(string(mooring_no) + '_detided');
ty = slp.datetime(ind);
y = slp.(string(mooring_no))(ind);

%%% Wavelet coherence
axes(ax4);
plot_wavelet_coherence(ax4, [], tx, x, y, 'PlotPhase', false, 'RunMonteCarlo', true, ...
    'NSurrogates', 200, 'SigLevel', 95, 'YTicks', [1 4 16 64 256 1024]);
hold on
xline(yearfrac(datenum([t1 t2])), 'r', 'LineWidth', lw)
xlabel(' ')
set(gca, 'FontSize', fs)
text(0.02, 0.93, '(d)', 'Units', 'normalized', 'Color', 'r',...
    'HorizontalAlignment', 'center', 'FontSize', fs, ...
    'FontWeight', 'bold', 'BackgroundColor', 'none', 'Parent', ax4)

