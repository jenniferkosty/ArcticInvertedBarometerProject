function [fave, coherence, cospectrum, transfer_function, phase, theor_coherence_crit, coherence_crit, cospectrum_crit, transfer_function_crit] = coherence_phase_pwelch(t, x, y, varargin)

p = inputParser;
p.FunctionName = mfilename;

addRequired(p,'t',@(v) isdatetime(v) && isvector(v));
addRequired(p,'x',@(v) isnumeric(v) && isvector(v) && numel(v)==numel(t));
addRequired(p,'y',@(v) isnumeric(v) && isvector(v) && numel(v)==numel(t));

addParameter(p,'NumSegments', 8,@(v) (isnumeric(v) && isvector(v)));
addParameter(p,'LineStyle','-',@(s) ischar(s) || isstring(s));
addParameter(p,'Color',rand(1,3),@(c) (ischar(c) || isstring(c)) || (isnumeric(c) && numel(c)==3));
addParameter(p,'Label', ' ',@(c) (ischar(c) || isstring(c)));
addParameter(p,'PlotPhase', 0,@(v) (isnumeric(v) && isvector(v)));
addParameter(p,'PlotCoherence', 0,@(v) (isnumeric(v) && isvector(v)));
addParameter(p,'PlotCrossSpectrum', 0,@(v) (isnumeric(v) && isvector(v)));
addParameter(p,'PlotTransferFunction', 0,@(v) (isnumeric(v) && isvector(v)));

% --- surrogate options
addParameter(p, 'MonteCarloSignificance', 0, @(v) (isnumeric(v) && isvector(v)))
addParameter(p,'Nsurr', 50, @(v) isnumeric(v) && isscalar(v) && v>=10);
addParameter(p,'Alpha', 0.05, @(v) isnumeric(v) && isscalar(v) && v>0 && v<1);

parse(p, t, x, y, varargin{:});
opt = p.Results;

x = x(:);
y = y(:);

%%% Finding near-continuous sections of data
good = isfinite(x) & isfinite(y);
starts = find(diff([0; good]) == 1);
ends   = find(diff([good; 0]) == -1);
[~, maxlength] = max(abs(starts - ends));
ind = starts(maxlength):ends(maxlength);
%ind = starts(1):ends(end);
x = x(ind);
y = y(ind);

%%% Removing mean
x = x - mean(x);
y = y - mean(y);

if length(x) < 20
    fave = NaN;
    coherence = NaN;
    cospectrum = NaN;
    transfer_function = NaN;
    phase = NaN;
    theor_coherence_crit = NaN;
    coherence_crit = NaN;
    cospectrum_crit = NaN;
    transfer_function_crit = NaN;
    return
end

%%% Computing estimate
dt = (datenum(t(2)) - datenum(t(1))); % in days
Fs = 1/dt;
segmentLength = fix(length(x)/opt.NumSegments); %%% Number of segments
N = length(x);
noverlap = floor(0.5*segmentLength); %%% Overlap between segments (typically 50-75%)
[Pxx, f] = pwelch(x,segmentLength,noverlap,N,Fs);
[Pyy, ~] = pwelch(y,segmentLength,noverlap,N,Fs);
[Pxy, ~] = cpsd(x, y,segmentLength,noverlap,N,Fs);

%%% Theoretical significance level for coherence based on ensemble number
K = 1+(N-segmentLength)/(segmentLength-noverlap);
alpha = opt.Alpha;
theor_coherence_crit = 1-alpha^(1/(K-1));

%%% Averaging in log-space
ff = f;
avgfrom = 1:length(ff);
a = linspace(1,log10(avgfrom(end)),100);
b = round(10.^a);
nuse = diff(b);
nix=[0 cumsum(nuse)];
nb   = length(nix)-1;

xxave = nan(1,length(nix)-1);
yyave = nan(1,length(nix)-1);
xyave = nan(1,length(nix)-1);
fave=nan*ones(1,length(nix)-1);
for n=1:nb
    ii = (nix(n)+1):nix(n+1);
    if isempty(ii)
        ii = nix(n)+1;
    end
    xxave(n) = mean(Pxx(ii), 'omitnan');
    yyave(n) = mean(Pyy(ii), 'omitnan');
    xyave(n) = mean(Pxy(ii), 'omitnan');
    fave(n)=   mean(ff(ii), 'omitnan');
end
[fave, idx] = unique(fave);
xxave = xxave(idx);
yyave = yyave(idx);
xyave = xyave(idx);

%%% Computing cross-spectrum, coherence, phase
phase = (angle(xyave));  
coherence = (abs(xyave).^2) ./ (xxave .* yyave);
cospectrum = (abs(xyave));

%%% Computing transfer function
transfer_function   = abs(xyave ./ xxave);

%%% Frequency-dependent significance testing using Monte-Carlo re-sampling
coherence_crit = NaN(size(coherence));
cospectrum_crit = NaN(size(cospectrum));
transfer_function_crit = NaN(size(transfer_function));
if opt.MonteCarloSignificance == 1

    Ns = opt.Nsurr;
    coherence_surr = nan(Ns, length(coherence));
    cospectrum_surr = nan(Ns, length(cospectrum));
    transfer_function_surr = nan(Ns, length(transfer_function));

    %%% Looping through surrogates
    for s = 1:Ns

        %%% Generating surrogate for y (while preserving power spectrum)
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

        %%% Computing psd, csd
        [PyyS, ~ ] = pwelch(yS, segmentLength, noverlap, N, Fs);
        [PxyS, ~ ] = cpsd(x, yS, segmentLength, noverlap, N, Fs);

        xyaveS = nan(1, nb);
        yyaveS = nan(1, nb);

        %%% Spectral averaging
        for n = 1:nb
            ii = (nix(n)+1):nix(n+1);
            if isempty(ii)
                ii = nix(n)+1;
            end
            yyaveS(n) = mean(PyyS(ii), 'omitnan');
            xyaveS(n) = mean(PxyS(ii), 'omitnan');
        end
        yyaveS = yyaveS(idx);
        xyaveS = xyaveS(idx);

        %%% Computing coherence
        coherence_surr(s,:) = (abs(xyaveS).^2) ./ (xxave .* yyaveS);

        %%% Computing cospectrum
        cospectrum_surr(s,:) = (abs(xyaveS));

        %%% Computing transfer function
        transfer_function_surr(s,:)   = abs(xyaveS ./ xxave);
    end

    %%% Computing significance level for each frequency
    q = 1 - opt.Alpha;
    coherence_crit = quantile(coherence_surr, q, 1);
    cospectrum_crit = quantile(cospectrum_surr, q, 1);
    transfer_function_crit = quantile(transfer_function_surr, q, 1);
end

lw = 2; % Linewidth

if opt.PlotPhase == 1

    %%% Plotting
    h = plot(fave,  phase, 'LineWidth', lw, 'Color', opt.Color, ...
        'LineStyle', opt.LineStyle);
    hold on

    %%% Adding label
    h(1).DisplayName = opt.Label;

    yticks(-180:60:180)
    xlog;
    xlabel('Frequency [cycles/day]')
    ylabel("Phase [deg]")
    grid on
    axis tight

    %%% Adding legend
    if ~strcmp(opt.Label, ' ')
        legend('show')
    end

elseif opt.PlotCoherence == 1

    %%% Plotting
    hold on
    
    h = plot(fave,  coherence, 'LineWidth', lw, 'Color', opt.Color, ...
        'LineStyle', opt.LineStyle);
    %plot(f, abs(Pxy).^2 ./ (Pxx .* Pyy), 'LineWidth', 0.5, 'Color', [0.5 0.5 0.5])
   
    %%% Adding Label
    h.DisplayName = opt.Label;

    yticks(0:0.1:1)
    xlog;
    xlabel('Frequency [cycles/day]')
    ylabel("Magnitude-Squared Coherence")
    grid on
    axis tight

    %%% Adding legend
    if ~strcmp(opt.Label, ' ')
        legend('show')
    end

elseif opt.PlotCrossSpectrum == 1

    %%% Plotting
    hold on
    
    h = plot(fave,  cospectrum, 'LineWidth', lw, 'Color', opt.Color, ...
        'LineStyle', opt.LineStyle);
   

    %%% Adding Label
    h(1).DisplayName = opt.Label;
    xlog;
    xlabel('Frequency [cycles/day]')
    ylabel("Cross-Spectrum")
    grid on
    axis tight

    %%% Adding legend
    if ~strcmp(opt.Label, ' ')
        legend('show')
    end

elseif opt.PlotTransferFunction == 1

    %%% Plotting
    hold on
    
    h = plot(fave, transfer_function, 'LineWidth', lw, 'Color', opt.Color, ...
        'LineStyle', opt.LineStyle);
   
    %%% Adding Label
    h(1).DisplayName = opt.Label;
    xlog;
    xlabel('Frequency [cycles/day]')
    ylabel("Transfer Function")
    grid on
    axis tight

    %%% Adding legend
    if ~strcmp(opt.Label, ' ')
        legend('show')
    end


end

end


