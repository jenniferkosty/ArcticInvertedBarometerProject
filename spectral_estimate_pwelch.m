function [fave, xave, P_theory, P_crit] = spectral_estimate_pwelch(t, x, varargin)

p = inputParser;
p.FunctionName = mfilename;

addRequired(p,'t',@(v) isdatetime(v) && isvector(v));
addRequired(p,'x',@(v) isnumeric(v) && isvector(v) && numel(v)==numel(t));

addParameter(p,'NumSegments', 8, @(v) isnumeric(v) && isscalar(v));
addParameter(p,'LineStyle', '-', @(s) ischar(s) || isstring(s));
addParameter(p,'Color', rand(1,3), @(c) (ischar(c) || isstring(c)) || (isnumeric(c) && numel(c)==3));
addParameter(p,'Label', ' ', @(c) ischar(c) || isstring(c));
addParameter(p,'PlotUnsmoothed', 0, @(v) isnumeric(v) && isscalar(v));
addParameter(p,'MakeFigure', 0, @(v) isnumeric(v) && isscalar(v));
addParameter(p,'PlotRedNoise', 0, @(v) isnumeric(v) && isscalar(v));
addParameter(p,'SignificanceLevel', 0.05, @(v) isnumeric(v) && isscalar(v) && v > 0 && v < 1);
addParameter(p,'ExcludeBand', [], @(v) isempty(v) || (isnumeric(v) && numel(v)==2));
addParameter(p,'FitFreqRange', [], @(v) isempty(v) || (isnumeric(v) && numel(v)==2));
addParameter(p,'BackgroundModel', 'powerlaw', @(s) any(strcmp(s, {'powerlaw','ar1'})));
addParameter(p,'AR1', [], @(v) isempty(v) || (isnumeric(v) && isscalar(v) && v > -1 && v < 1));

parse(p, t, x, varargin{:});
opt = p.Results;

%%% Formatting data
x = x(:);
t = t(:);

%%% Finding near-continuous sections of data
good = isfinite(x);
starts = find(diff([0; good]) == 1);
ends   = find(diff([good; 0]) == -1);
[~, maxlength] = max(abs(starts - ends));
ind = starts(maxlength):ends(maxlength);
x = x(ind);
if length(x) < 20
    fave     = NaN;
    xave     = NaN;
    P_theory = NaN;
    P_crit   = NaN;
    return
end

%%% Removing mean
x = x - mean(x, 'omitnan');

%%% Computing pwelch estimate
dt           = datenum(t(2)) - datenum(t(1));   % days
Fs           = 1/dt;                             % cycles/day
segmentLength = fix(length(x)/opt.NumSegments);
N            = length(x);
noverlap     = floor(segmentLength/2);
[Pxx, freq]  = pwelch(x, segmentLength, noverlap, N, Fs);

% =========================================================================
% BACKGROUND MODEL
% =========================================================================

ff       = freq(2:end);      % all non-DC frequencies
Pxx_nodc = Pxx(2:end);

if strcmp(opt.BackgroundModel, 'powerlaw')

    % -------------------------------------------------------------------
    % POWER LAW BACKGROUND:  P(f) = A * f^(-beta)
    %
    % Fit a straight line in log-log space:
    %   log10(P) = log10(A) - beta * log10(f)
    %
    % The fit can be restricted to a frequency range where the power law
    % is a good description of the background (FitFreqRange), and a
    % suspected peak band can be masked out (ExcludeBand).
    % -------------------------------------------------------------------

    log_f = log10(ff);
    log_P = log10(Pxx_nodc);

    fit_mask = isfinite(log_f) & isfinite(log_P);

    if ~isempty(opt.ExcludeBand)
        fit_mask = fit_mask & ...
            ~(ff >= opt.ExcludeBand(1) & ff <= opt.ExcludeBand(2));
    end

    if ~isempty(opt.FitFreqRange)
        fit_mask = fit_mask & ...
            ff >= opt.FitFreqRange(1) & ff <= opt.FitFreqRange(2);
    end

    coeffs        = polyfit(log_f(fit_mask), log_P(fit_mask), 1);
    beta          = -coeffs(1);
    P_theory_full = 10.^polyval(coeffs, log_f);

    background_label = sprintf('Power law (\\beta = %.2f)', beta);

elseif strcmp(opt.BackgroundModel, 'ar1')

    % -------------------------------------------------------------------
    % AR1 BACKGROUND:
    %
    %   P(f) = P0 * (1 - r1^2) / (1 - 2*r1*cos(2*pi*f*dt) + r1^2)
    %
    % where:
    %   r1  = lag-1 autocorrelation of the time series
    %   P0  = mean spectral power (scaling factor so the background
    %         integrates to the same total variance as the data)
    %   dt  = sampling interval in days
    %
    % If the user supplies 'AR1', that value is used directly.
    % Otherwise r1 is estimated from the data using the lag-1
    % autocorrelation of the observed time series.
    %
    % The AR1 spectrum has a maximum slope of f^-2 at high frequencies.
    % If your observed spectrum is steeper than this, consider using
    % the power law background instead.
    % -------------------------------------------------------------------

    if ~isempty(opt.AR1)
        % User-supplied r1
        r1 = opt.AR1;
    else
        % Estimate r1 from the data.
        %
        % The lag-1 autocorrelation is:
        %   r1 = sum(x(1:N-1) .* x(2:N)) / sum(x.^2)
        %
        % We demean first so r1 reflects temporal memory, not the mean.
        xd = x - mean(x, 'omitnan');
        r1 = sum(xd(1:end-1) .* xd(2:end)) / sum(xd.^2);
        r1 = max(-0.999, min(0.999, r1));   % clamp to valid range
    end

    % Theoretical AR1 spectral shape (unnormalised)
    %
    %   The denominator of the AR1 spectrum:
    %       D(f) = 1 - 2*r1*cos(2*pi*f*dt) + r1^2
    %
    %   The (1 - r1^2) numerator ensures the spectrum integrates to
    %   unit variance for a unit-variance process. We then scale by
    %   P0 so the background matches the observed spectral level.
    %
    %   P0 is chosen so that the mean of P_theory_full equals the
    %   mean of the observed spectrum over the fit range (or over
    %   all frequencies if no FitFreqRange is given). This is a
    %   least-squares scaling in linear space.

    D             = 1 - 2*r1*cos(2*pi*ff*dt) + r1^2;
    ar1_shape     = (1 - r1^2) ./ D;    % theoretical shape, unit variance

    % Build scaling mask (same logic as power law fit mask)
    scale_mask = isfinite(Pxx_nodc) & isfinite(ar1_shape);

    if ~isempty(opt.ExcludeBand)
        scale_mask = scale_mask & ...
            ~(ff >= opt.ExcludeBand(1) & ff <= opt.ExcludeBand(2));
    end

    if ~isempty(opt.FitFreqRange)
        scale_mask = scale_mask & ...
            ff >= opt.FitFreqRange(1) & ff <= opt.FitFreqRange(2);
    end

    % Least-squares scalar to match observed spectral level:
    %   P0 = sum(Pxx * ar1_shape) / sum(ar1_shape^2)
    P0            = sum(Pxx_nodc(scale_mask) .* ar1_shape(scale_mask)) / ...
                    sum(ar1_shape(scale_mask).^2);
    P_theory_full = P0 * ar1_shape;

    background_label = sprintf('AR1 (r_1 = %.2f)', r1);

end

% =========================================================================
% CHI-SQUARED SIGNIFICANCE TEST
% =========================================================================
%
% The pwelch estimate at each frequency follows a scaled chi-squared
% distribution with dof = 2K degrees of freedom, where K is the number
% of segments contributing to the average.
%
% The (1-alpha) significance threshold is:
%   P_crit(f) = P_theory(f) * chi2inv(1-alpha, dof) / dof

K             = 1 + (N - segmentLength) / (segmentLength - noverlap);
dof           = 2 * K;
alpha         = opt.SignificanceLevel;
P_crit_full   = P_theory_full * chi2inv(1 - alpha, dof) / dof;

% =========================================================================
% LOG-SPACE AVERAGING
% =========================================================================

avgfrom = 1:length(ff);
a       = linspace(1, log10(avgfrom(end)), 100);
b       = round(10.^a);
nuse    = diff(b);
nix     = [0 cumsum(nuse)];
nb      = length(nix) - 1;

xave         = nan(1, nb);
fave         = nan(1, nb);
P_theory_ave = nan(1, nb);
P_crit_ave   = nan(1, nb);

for n = 1:nb
    ii = (nix(n)+1):nix(n+1);
    if isempty(ii)
        ii = nix(n)+1;
    end
    xave(n)         = mean(Pxx_nodc(ii),     'omitnan');
    fave(n)         = mean(ff(ii),            'omitnan');
    P_theory_ave(n) = mean(P_theory_full(ii), 'omitnan');
    P_crit_ave(n)   = mean(P_crit_full(ii),   'omitnan');
end

[fave, idx] = unique(fave);
xave        = xave(idx);
P_theory    = P_theory_ave(idx);
P_crit      = P_crit_ave(idx);

% =========================================================================
% PLOTTING
% =========================================================================

if opt.MakeFigure == 1
    hold on;
    lw = 2;

    h = plot(fave, xave, 'LineWidth', lw, 'Color', opt.Color, ...
        'LineStyle', opt.LineStyle);
    h.DisplayName = opt.Label;
    if strcmp(opt.Label, ' ')
        h.HandleVisibility = 'off';
    end

    %%% Plotting unsmoothed version
    if opt.PlotUnsmoothed == 1
        pu = plot(ff, Pxx_nodc, 'LineWidth', 1, 'Color', [0.5 0.5 0.5], ...
            'LineStyle', '-');
        pu.HandleVisibility = 'off';
    end

    %%% Plotting red noise
    if opt.PlotRedNoise == 1
        h_theory = plot(fave, P_theory, 'k--', 'LineWidth', 1.5);
        h_theory.DisplayName = background_label;

        h_crit = plot(fave, P_crit, 'r--', 'LineWidth', 1.5);
        h_crit.DisplayName = sprintf('%d%% significance', round((1-alpha)*100));
    end

    ylabel('Power Spectral Density')
    xlabel('Frequency [cycles/day]')
    set(gca, 'Xscale', 'log');
    set(gca, 'Yscale', 'log');
    grid on

    %%% Adding legend if appropriate
    if ~strcmp(opt.Label, ' ') || opt.PlotRedNoise == 1
        legend('show');
    end
end
end