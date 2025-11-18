clear all; clc; close all;

% Wave parameters
g       = 9.81;
h       = Inf;
Hs      = 1.5; 
T0      = 8;
w0      = 2*pi/T0;
gamma   = 3.3;

% Directional spreading
beta = deg2rad(45);
s_theta = 5;

% Simulation
Nw = 50;
Nth = 50;
w_min = 0.2*w0;
w_max = 3*w0;
N = 150;
w = linspace(w_min,w_max, N);
dw = w(2) - w(1);

theta = linspace(beta - pi/2, beta + pi/2, Nth);
dth = theta(2) - theta(1);

% JONSWAP spectrum
S = wavespec(7,[Hs, w0, gamma], w, 1);

% Sea Surface elevation
T_end = 150;
dt = 0.1;
t = 0:dt:T_end;

dw = w(2)-w(1);
phi  = 2*pi*rand(size(w));
Ak = sqrt(2*S*dw);
zeta = sum(Ak(:).*cos(w(:).*t + phi(:)), 1);

% Directional spreading
D = cos((theta - beta)/2).^(2*s_theta);
D(abs(theta-beta) > pi) = 0;
D = D / (sum(D)*dth);

% Grid
Lx = 300;  Ly = 300;
Nx = 160;  Ny = 160;
x  = linspace(0, Lx, Nx);
y  = linspace(0, Ly, Ny);
[X, Y] = meshgrid(x, y);

% -------------------- Component wavenumbers & phases ---------------------
k = w2k(w, h, g);                      % solve dispersion (deep or finite depth)
rng(1);                                 % reproducible
phi = 2*pi*rand(Nw, Nth);               % random phases

% Component amplitudes from 2D spectrum: A = sqrt(2 S(w) D(th) dw dth)
A = sqrt(2 * (S(:) * D(:).') * dw * dth);   % Nw x Nth

% -------------------- Sea-surface elevation snapshot ----------------------
t0   = 50;
eta = zeros(Ny, Nx);
for i = 1:Nw
    for j = 1:Nth
        kx  = k(i)*cos(theta(j));
        ky  = k(i)*sin(theta(j));
        eta = eta + A(i,j) * cos(kx*X + ky*Y - w(i)*t0 + phi(i,j));
    end
end

%% Plot
% figure('Color','w');
% plot(w, S, 'LineWidth', 1.6); grid on;
% xline(w0,'--','\omega_0'); xlabel('\omega [rad/s]');
% ylabel('S_\eta(\omega) [m^2 s]'); title('JONSWAP spectrum (normalized to H_s)');

figure;
plot(t, zeta, 'LineWidth', 1);
grid on;
xlabel('t [s]'); ylabel('\zeta [m]');
title('Realization of long-crested sea surface elevation');

figure('Color','w');
surf(x, y, eta, 'EdgeColor', [0.3 0.3 0.3], 'EdgeAlpha', 0.35);
shading interp; colormap jet; view(-32,30); grid on;
xlabel('x [m]'); ylabel('y [m]'); zlabel('\eta [m]');
title(sprintf('Short-crested sea realization, %d components', Nw*Nth));


function k = w2k(w, h, g)
% Solve dispersion relation: w^2 = g k tanh(k h)  (deep water: tanh -> 1)
    k = w.^2./g;                        % deep-water initial guess
    if isfinite(h)
        for iter = 1:30                 % Newton–Raphson
            f  = g.*k.*tanh(k*h) - w.^2;
            df = g.*(tanh(k*h) + k.*h.*sech(k*h).^2);
            k  = k - f./df;
        end
    end
end