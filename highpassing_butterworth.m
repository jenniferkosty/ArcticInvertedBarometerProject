function [hp, lp] = highpassing_butterworth(x, sampling_interval, cutoff)

% INPUT:
% x, time series
% sampling_interval, in minutes
% cutoff, in hours

x = x(:);

dt = 60*sampling_interval;          % seconds per sample
fs = 1/dt;                          % sampling frequency in Hz

cutoff_period = cutoff * 3600;      % convert to seconds
cutoff_freq = 1 / cutoff_period;    % Hz

order = 4;
Wn = cutoff_freq / (fs/2);          % normalize to Nyquist

[b, a] = butter(order, Wn, 'high');
x = x(:);
hp = NaN(size(x));

% Find continuous blocks of good data
good = isfinite(x);
starts = find(diff([0; good]) == 1);
ends   = find(diff([good; 0]) == -1);
segment_lengths = ends - starts;
min_length = order * 10;  % need enough points to filter

% High-passing individual segments
for i = 1:length(starts)
    idx = starts(i):ends(i);
    if length(idx) > min_length
        hp(idx) = filtfilt(b, a, x(idx));
    end
end

lp = x - hp;

make_figure = 0;
if make_figure == 1

    x = x(:);
    hp = hp(:);
    lp = lp(:);

    figure('Position', [10 10 1000 500])
    fs = 12; lw = 2;
    tiledlayout(1,2)
    sgtitle('Butterworth, Order = ' + string(order))

    %%% Plotting original and highpassed time series
    nexttile()
    plot(x)
    hold on
    plot(hp)
    xlabel('Time')
    grid on
    set(gca, 'FontSize', fs)

    %%% Plotting spectral estimates
    nexttile()
    [~, ii] = max(segment_lengths);
    dt = sampling_interval/60;
    [psio,lambda]=sleptap(length(x(starts(ii):ends(ii))),16);
    [f,s]=mspec(dt,x(starts(ii):ends(ii)),psio,'cyclic'); %Spectrum of original
    [~,sr]=mspec(dt,hp(starts(ii):ends(ii)),psio,'cyclic');%Spectrum of residual

    % [psio,lambda]=sleptap(length(x),16);
    % x(isnan(x)) = 0;
    % hp(isnan(hp)) = 0;
    % [f,s]=mspec(dt,x,psio,'cyclic'); %Spectrum of original
    % [~,sr]=mspec(dt,hp,psio,'cyclic');%Spectrum of residual

    a(1) = plot(f,s, 'LineWidth', lw, 'DisplayName', 'Original');
    hold on
    a(2) = plot(f, sr, 'LineWidth', lw, 'DisplayName', 'Highpassed');
    xlabel('Frequency (cycles/hour)')
    ylog, axis tight,ax=axis;
    grid on
    set(gca, 'FontSize', fs)
    legend()

end

end