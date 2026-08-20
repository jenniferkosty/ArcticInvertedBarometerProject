function [f, coherence, cospectrum, transfer_function, phase, theor_coherence_crit, coherence_crit, cospectrum_crit, transfer_function_crit] = coherence_phase_pwelch_full_record(tx_full, x_full, ty_full, y_full, varargin)

p = inputParser;
p.FunctionName = mfilename;

addRequired(p,'tx_full',@(v) isdatetime(v) && isvector(v));
addRequired(p,'x_full',@(v) isnumeric(v) && isvector(v) && numel(v)==numel(tx_full));
addRequired(p,'ty_full',@(v) isdatetime(v) && isvector(v));
addRequired(p,'y_full',@(v) isnumeric(v) && isvector(v) && numel(v)==numel(ty_full));

% optional parameters
addParameter(p, 'StartDates', tx_full(1), @(v) isdatetime(v) && isvector(v));
addParameter(p, 'EndDates', tx_full(end), @(v) isdatetime(v) && isvector(v));
addParameter(p,'NumSegments', 8,@(v) (isnumeric(v) && isvector(v)));
addParameter(p, 'MonteCarloSignificance', 0, @(v) (isnumeric(v) && isvector(v)))
addParameter(p,'Nsurr', 50, @(v) isnumeric(v) && isscalar(v) && v>=10);
addParameter(p,'Alpha', 0.05, @(v) isnumeric(v) && isscalar(v) && v>0 && v<1);

parse(p, tx_full, x_full, ty_full, y_full, varargin{:});
opt = p.Results;

x_full = x_full(:);
y_full = y_full(:);

dt = (datenum(tx_full(2)) - datenum(tx_full(1))); % in days
Fs = 1/dt;
nYears = length(opt.StartDates);

%%% Creating frequency grid for averaging
a = 1/365;
if dt == hours(6)
    b = 2;
elseif dt == days(1)
    b = 1/2;
end
Nf = 100;
f = logspace(log10(a), log10(b), Nf);

%%% Looping through annual segments
Pxx_interp = NaN(length(opt.StartDates), length(f));
Pyy_interp = NaN(length(opt.StartDates), length(f));
Pxy_interp = complex(NaN(length(opt.StartDates), length(f)));
x_segments = cell(1, nYears);
y_segments = cell(1, nYears);
K_total = 0;

for i = 1:nYears

    %%% Isolating data for a given segment
    idx = tx_full >= opt.StartDates(i) & tx_full <= opt.EndDates(i);
    x = x_full(idx);
    idy = ty_full >= opt.StartDates(i) & ty_full <= opt.EndDates(i);
    y = y_full(idy);
    if length(x) ~= length(y)
        continue
    end

    %%% Finding near-continuous sections of data
    good = isfinite(x) & isfinite(y);
    starts = find(diff([0; good]) == 1);
    ends   = find(diff([good; 0]) == -1);
    [~, maxlength] = max(abs(starts - ends));
    ind = starts(maxlength):ends(maxlength);
    x = x(ind);
    y = y(ind);

    %%% Removing mean
    x = x - mean(x, 'omitnan');
    y = y - mean(y, 'omitnan');

    if length(x) < 20
        continue
    end


    segmentLength = fix(length(x)/opt.NumSegments); %%% Number of segments
    N = length(x);
    noverlap = floor(0.5*segmentLength); %%% Overlap between segments (typically 50-75%)

    %%% Store for surrogate loop
    x_segments{i}          = x;
    y_segments{i}          = y;
    seg_params(i).N             = N;
    seg_params(i).segmentLength = segmentLength;
    seg_params(i).noverlap      = noverlap;

    %%% Computing estimates
    [Pxx, ff] = pwelch(x,segmentLength,noverlap,N,Fs);
    [Pyy, ~] = pwelch(y,segmentLength,noverlap,N,Fs);
    [Pxy, ~] = cpsd(x, y,segmentLength,noverlap,N,Fs);

    %%% Interpolating to common frequency grid
    Pxx_interp(i,:) = interp1(ff, Pxx, f);
    Pyy_interp(i,:) = interp1(ff, Pyy, f);
    Pxy_interp(i,:) = interp1(ff, real(Pxy), f) + ...
        1i* interp1(ff, imag(Pxy), f);

    %%% Accumulating total degrees of freedom across all deployments
    K_i = 1+(N-segmentLength)/(segmentLength-noverlap);
    K_total = K_total + K_i;
end

%%% Theoretical significance level
alpha = opt.Alpha;
theor_coherence_crit = 1-alpha^(1/(K_total-1));

%%% Averaging across deployments
Pxx_mean = mean(Pxx_interp, 1, 'omitnan');
Pyy_mean = mean(Pyy_interp, 1, 'omitnan');
Pxy_mean = mean(Pxy_interp, 1, 'omitnan');

%%% Computing cross-spectrum, coherence, phase
phase = (angle(Pxy_mean));  
coherence = (abs(Pxy_mean).^2) ./ (Pxx_mean .* Pyy_mean);
cospectrum = (abs(Pxy_mean));

%%% Computing transfer function
transfer_function   = abs(Pxy_mean) ./ Pxx_mean;

%%% Frequency-dependent significance testing using Monte-Carlo re-sampling
coherence_crit = NaN(1, Nf);
cospectrum_crit = NaN(1, Nf);
transfer_function_crit = NaN(1, Nf);

if opt.MonteCarloSignificance == 1

    Ns = opt.Nsurr;
    coherence_surr = nan(Ns, Nf);
    cospectrum_surr = nan(Ns, Nf);
    transfer_function_surr = nan(Ns, Nf);

    %%% Looping through surrogates
    for s = 1:Ns

        PxxS_interp = NaN(nYears, Nf);
        PyyS_interp = NaN(nYears, Nf);
        PxyS_interp = complex(NaN(nYears, Nf));

        for i = 1:nYears

            x = x_segments{i};
            y = y_segments{i};
            if length(x) < 20 % year was skipped above
                continue; 
            end    

            N             = seg_params(i).N;
            segmentLength = seg_params(i).segmentLength;
            noverlap      = seg_params(i).noverlap;

            %%% Phase-randomize y
            Y           = fft(y);
            Nhalf       = floor(N/2);
            rand_phases = exp(sqrt(-1) * 2 * pi * rand(N, 1));

            rand_phases(1) = 1;                                            % preserve DC

            if mod(N, 2) == 0
                rand_phases(Nhalf+1) = real(rand_phases(Nhalf+1));   % Nyquist (even N only)
            end
            mirror_idx = Nhalf+2 : N;

            if mod(N,2) == 0
                source_idx = Nhalf : -1 : 2;
            else
                source_idx = Nhalf+1 : -1 : 2;
            end
            rand_phases(mirror_idx) = conj(rand_phases(source_idx));       % conjugate symmetry
            yS = real(ifft(Y .* rand_phases));

            %%% Surrogate spectra
            [PxxS, ff ] = pwelch(x,  segmentLength, noverlap, N, Fs);
            [PyyS, ~ ] = pwelch(yS,  segmentLength, noverlap, N, Fs);
            [PxyS, ~ ] = cpsd(x, yS, segmentLength, noverlap, N, Fs);

            %%% Interpolate
            PxxS_interp(i,:) =       interp1(ff, PxxS,       f, 'linear', NaN);
            PyyS_interp(i,:) =       interp1(ff, PyyS,       f, 'linear', NaN);
            PxyS_interp(i,:) = interp1(ff, real(PxyS), f, 'linear', NaN) + ...
                           1i* interp1(ff, imag(PxyS), f, 'linear', NaN);
        end

        %%% Average surrogate spectra across years
        PxxS_mean = mean(PxxS_interp, 1, 'omitnan');
        PyyS_mean = mean(PyyS_interp, 1, 'omitnan');
        PxyS_mean = mean(PxyS_interp, 1, 'omitnan');

        %%% Surrogate derived quantities
        coherence_surr(s,:)         = abs(PxyS_mean).^2 ./ (PxxS_mean .* PyyS_mean);
        cospectrum_surr(s,:)        = abs(PxyS_mean);
        transfer_function_surr(s,:) = abs(PxyS_mean) ./ PxxS_mean;

    end

    %%% Significance thresholds
    q                      = 1 - opt.Alpha;
    coherence_crit         = quantile(coherence_surr,         q, 1);
    cospectrum_crit        = quantile(cospectrum_surr,        q, 1);
    transfer_function_crit = quantile(transfer_function_surr, q, 1);

end

end



