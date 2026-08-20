function harmonics = computing_harmonics(cyc, x, t)
   
%%% INPUT:
% cyc: array with the desired harmonics, in terms of # of cycles per year
% x: time series
% t: time array in datetime

%%% Formatting harmonic cycles array
cyc = cyc(:);
harmonics.cycles = cyc;

%%% Formatting time array
t = t(:);

%%% Formatting data array
x = x(:);

%%% Removing NaNs
nan_mask = ~isnan(x);
x_orig = x;          % keep original for reinsertion
t_orig = t;          % keep original time array
x = x(nan_mask);
t = t(nan_mask);

%%% Converting time to fractional years
t_years = years(t - t_orig(1));

%%% Building design matrix
n = length(t);
m = length(cyc);
X = ones(n, 1 + 2*m);
col = 2;
for i = 1:m
    k = cyc(i);
    X(:,col)   = cos(2*pi*k*t_years);
    X(:,col+1) = sin(2*pi*k*t_years);
    col = col + 2;
end

%%% Least squares fit
beta = X \ x;

%%% Reconstruct signal
x_fit = X * beta;

%%% Calculating variance explained
residual = x - x_fit;
harmonics.total_var_explained = 100*(1 - var(residual,'omitnan')/var(x,'omitnan')); % percent

%%% Calculating variance explained per harmonic
var_total = var(x);
var_full  = var(residual);
VE = zeros(m,1);
for i = 1:m
    
    %%% columns for this harmonic in X
    cols = [2*i, 2*i+1];
    
    %%% design matrix WITHOUT harmonic i
    X_red = X;
    X_red(:,cols) = [];
    
    %%% fit reduced model
    beta_red = X_red \ x;
    res_red = x - X_red*beta_red;
    
    %%% incremental variance explained
    VE(i) = 100 * (var(res_red) - var_full) / var_total;
end
harmonics.var_each = VE;

%%% Extracting amplitude and phases
a   = zeros(m,1);
b   = zeros(m,1);
A   = zeros(m,1);
phi = zeros(m,1);
idx = 2;
for i = 1:m
    a(i) = beta(idx);
    b(i) = beta(idx+1);
    A(i) = hypot(a(i), b(i));
    phi(i) = atan2(b(i), a(i));
    idx = idx + 2;
end

%%% Saving data
harmonics.a = a;
harmonics.b = b;
harmonics.A = A;
harmonics.phi = phi;

%%% Reconstructing each harmonic for future plotting
n_orig = length(x_orig);
components_full = NaN(n_orig, m);
for i = 1:m
    k = cyc(i);
    comp = a(i)*cos(2*pi*k*t_years) + b(i)*sin(2*pi*k*t_years);
    components_full(nan_mask, i) = comp;
end

%%% Saving data
x_fit_full = NaN(n_orig, 1);
x_fit_full(nan_mask) = x_fit;

harmonics.components   = components_full;
harmonics.time         = t_orig;
harmonics.reconstruction = x_fit_full;
harmonics.mean         = beta(1);

end

