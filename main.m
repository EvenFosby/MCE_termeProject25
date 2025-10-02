 % --- User inputs ---
g     = 9.81;
Hs    = 4.0;       % [m]
Tp    = 10.0;      % [s]
gamma = 3.3;       % JONSWAP peak enhancement
Tsim  = 300;      % [s]
dt    = 0.1;       % [s]
deep  = true;      % set false if you want finite depth (then provide h)
h     = 100.0;     % [m], only used if deep=false

% --- Grids ---
t  = (0:dt:Tsim).';
dw = 2*pi/Tsim;
wmax = min(5*(2*pi/Tp), pi/dt*0.9);   % keep under Nyquist
w  = (1:floor(wmax/dw)).' * dw;       % start at k=1 (no DC); size Nw x 1
Nw = numel(w);

% --- JONSWAP shape (unscaled) ---
wp = 2*pi/Tp;
sigma = 0.09*ones(Nw,1); sigma(w<=wp) = 0.07;
Phi = (g^2)./w.^5 .* exp(-1.25*(wp./w).^4) .* gamma.^( exp( - (w-wp).^2 ./ (2*(sigma.^2)*(wp^2)) ) );

% --- Scale to match Hs ---
m0_shape = trapz(w, Phi);
alpha = (Hs^2/16) / m0_shape;
S = alpha * Phi;          % S_eta(w) [m^2 s]

% --- Synthesis: 1-point time series ---
rng(1);                   % for reproducibility
phi = 2*pi*rand(Nw,1);
Ak  = sqrt(2*S*dw);       % amplitudes
Wt  = t*w.';              % (Nt x Nw)
eta = cos(Wt + phi.');    % (Nt x Nw)
eta = eta * Ak;           % (Nt x 1) linear combination

% --- QC ---
eta_std = std(eta);
fprintf('Target std = %.3f m, realized std = %.3f m\n', Hs/sqrt(16), eta_std);

% --- OPTIONAL: directional, spatial field ---
useDirectional = false;  % set true to build eta(x,y,t)
if useDirectional
    % Direction grid + spreading
    Nb = 73;
    beta0 = deg2rad(0);  % mean direction coming from 0°
    beta = linspace(-pi, pi, Nb); db = beta(2)-beta(1);
    s = 10;              % spreading exponent (cos^2s)
    D_raw = cos((beta - beta0)/2).^(2*s); D_raw(D_raw<0) = 0;
    D = D_raw / trapz(beta, D_raw);  % normalize

    % Wavenumbers
    if deep
        k = (w.^2)/g;
    else
        % simple fixed-point iteration for dispersion: w^2 = g k tanh(k h)
        k = (w.^2)/g; % init deep
        for it=1:20
            k = (w.^2)./(g*tanh(k*h));
        end
    end

    % Spatial grid (small demo)
    x = linspace(-50,50,101); y = 0;  % 1D cut; extend to meshgrid for 2D
    [XX,TT,BB] = ndgrid(x, t, beta);  %#ok<ASGLU>
    % Component amplitudes
    Akn = sqrt(2*S*dw) .* sqrt(D*db); % (Nw x Nb), outer product below
    eta_xy_t = zeros(numel(t), numel(x));
    % Sum efficiently
    for iw = 1:Nw
        % phase term k*(x*cosb) - w*t + phi
        phase = (x.*cos(beta))'*k(iw) - w(iw)*t;
        % random phases per (w,beta)
        phin = 2*pi*rand(1,Nb);
        eta_xy_t = eta_xy_t + cos(phase + phin).* (Akn(iw,:));
    end
    % Example: take a time slice or point trace
    eta_x_t = eta_xy_t; % (Nt x Nx)
end

% After you build t, eta, w, S from the JONSWAP synthesis:
plot_wave_timeseries(t, eta, Hs);
plot_wave_spectrum(t, eta, w, S, 'title', 'JONSWAP check');

% If you form a directional spectrum S_wb = S(:).*D(:)' with beta:
plot_directional_spectrum(w, beta, S_wb, 'beta0', 0);

function plot_wave_timeseries(t, eta, Hs, varargin)
% PLOT_WAVE_TIMESERIES  Time trace with ±2σ lines and Hs markers.
% t [s], eta [m], Hs [m] (target).  Optional: 'title', 'ylim'
p = inputParser;
addParameter(p,'title','Wave elevation');
addParameter(p,'ylim',[]);
parse(p,varargin{:});
dt = mean(diff(t));
sig = std(eta);
Hs_est = 4*sig;

figure; 
plot(t, eta, 'LineWidth', 1); grid on;
hold on;
yline( 2*sig, '--', '\pm2\sigma'); 
yline(-2*sig, '--', ''); 
yline( Hs/4, ':', 'Hs/4 (target)', 'LabelVerticalAlignment','bottom');
yline(-Hs/4, ':', '', 'LabelVerticalAlignment','top');
xlabel('Time [s]'); ylabel('\eta [m]');
title(p.Results.title);

if ~isempty(p.Results.ylim), ylim(p.Results.ylim); end

txt = sprintf('dt = %.3g s | std = %.3g m | Hs(est)=%.3g m | Hs(target)=%.3g m', ...
               dt, sig, Hs_est, Hs);
annotation(gcf,'textbox',[.15 .8 .3 .12],'String',txt,'FitBoxToText','on','EdgeColor',[0.8 0.8 0.8]);
end

% -------------------------------------------------------------------------

function plot_wave_spectrum(t, eta, w_theory, S_theory, varargin)
% Overlay estimated PSD (Welch) vs theoretical S(ω).
p = inputParser;
addParameter(p,'title','Spectrum (estimated vs theory)');
addParameter(p,'nperseg',[]);
addParameter(p,'overlap',[]);
parse(p,varargin{:});

dt = mean(diff(t)); fs = 1/dt; N = numel(eta);
nper = p.Results.nperseg; if isempty(nper), nper = max(256, 2^floor(log2(N/8))); end
nover = p.Results.overlap; if isempty(nover), nover = floor(nper/2); end
step = nper - nover;
idx0 = 1:step:(N - nper + 1);
K = numel(idx0);
if K < 1, error('nperseg/overlap too large for the data length.'); end

win = hann_local(nper);
S_accum = 0;
for k = 1:K
    seg  = eta(idx0(k):idx0(k)+nper-1);
    segw = seg .* win;
    X = fft(segw, nper);
    % one-sided PSD in Hz
    P = (abs(X).^2) / (fs * sum(win.^2));
    P1 = P(1:floor(nper/2)+1);
    % double non-DC/Nyquist to preserve total power
    if nper > 2
        P1(2:end-1) = 2*P1(2:end-1);
    end
    S_accum = S_accum + P1;
end
Pxx = S_accum / K;                 % [m^2/Hz]
f    = (0:floor(nper/2))' * fs/nper;
w_est = 2*pi*f;                    % [rad/s]
S_est = Pxx / (2*pi);              % [m^2 s]

% QC
m0_est    = trapz(w_est, S_est);
m0_theory = trapz(w_theory, S_theory);
var_samp  = var(eta);

figure;
loglog(w_est, S_est, 'LineWidth', 1.2); hold on; grid on;
loglog(w_theory, S_theory, '--', 'LineWidth', 1.2);
xlabel('\omega [rad/s]'); ylabel('S_\eta(\omega) [m^2 s]');
legend('Estimated (Welch)', 'Theoretical', 'Location','best');
title(p.Results.title);

txt = sprintf('∫S_{est} d\\omega = %.3g (sample var=%.3g)\\n∫S_{theory} d\\omega = %.3g', ...
              m0_est, var_samp, m0_theory);
annotation(gcf,'textbox',[.14 .72 .34 .14],'String',txt,'FitBoxToText','on','EdgeColor',[0.8 0.8 0.8]);
end

function w = hann_local(N)
% Hanning window (toolbox-free)
w = 0.5 * (1 - cos(2*pi*(0:N-1)'/(N-1)));
end
% -------------------------------------------------------------------------

function plot_directional_spectrum(w, beta_grid, S_wb, varargin)
% PLOT_DIRECTIONAL_SPECTRUM  Heatmap of S(ω,β).
% w: [Nw x 1] rad/s, beta_grid: [1xNb] or [Nb x 1] rad, S_wb: [Nw x Nb]
p = inputParser;
addParameter(p,'beta0',NaN);
parse(p,varargin{:});

if size(S_wb,1) ~= numel(w)
    error('S_wb must be Nw x Nb to match w.');
end
beta_deg = beta_grid(:)' * 180/pi;

figure;
imagesc(beta_deg, w, S_wb); axis xy;
xlabel('\beta [deg]'); ylabel('\omega [rad/s]');
title('Directional spectrum S_\eta(\omega,\beta)');
cb = colorbar; cb.Label.String = '[m^2 s / rad]';
grid on;

if ~isnan(p.Results.beta0)
    hold on; xline(p.Results.beta0, 'w--', 'Mean dir');
end
end