close all; clear; clc; format default;

p = radarParams();
A = echoMatrix(p); % NxM samples


rawMatrix(A, p);
fprintf('dR = %.2f m | dv = %.2f m/s | Rmax = %.0f m | vmax = %.0f m/s | fD = %.0f Hz\n', ...
        p.dR, p.dv, p.Rmax, p.vmax, p.fD);

%% Parameters
% Initialize all parameters for simulation that we will need
function p = radarParams()
    p.c = 3*10^8; p.fc = 10e9;
    p.lambda = p.c/p.fc;
    p.B = 30e6;
    p.fs = 60e6;
    p.Tp = 5e-6; % Pulse width  (s), howlong transmitter on for each pulse
    p.PRF = 20e3; % pulses per second (1/PRI)
    p.PRI = 1/p.PRF; % pulse period (s)
    p.M = 64; % number of pulses we process togetherr for velocity
    % we set this to be high enough to get fine resolustion but not computationally too difficult 

    p.R0 = 2e3; % range of target at t = 0 (m)
    p.v = 30; % velocity of target (radial, constant positive means approaching) (m/s)
    p.A = 1; % echo amplitude at baseband, for noise measurements SNR
    
    p.N = round(p.PRI*p.fs); % number of samples per pulse interval (listening window: 3000)
    p.k = p.B/p.Tp; % chirp rate for LFM (1/s^2)
    p.fD = 2*p.v/p.lambda; % expected doppler shift from moving target (Hz), used to check results (using emag eqns)
    p.dR = p.c/(2*p.B); % range resolution of our radar set from bandwidth (m)
    p.dv = p.lambda/(2*p.M*p.PRI); % vel res (same as fro range res) (m/s)
    
    p.Rmax = p.c*p.PRI/2; % max range comes from PRI, past this range we get folding onto next Pulse(m)
    % min range is from Tp, we cant listen while tx is on, so something too
    % close will never be detected
    p.vmax = p.lambda*p.PRF/4; % max veel based on PRF and wavelength
    % past this aliasing occurs due to phase shift being too high(m/s)
    
    % Notice a trade off for Rmax and vmax for Rmax*vmax = c*lambda/8 = 1.125e6
    % increasing PRF decreases range max but increases vmax and vice versa.
    % chose values for PRI/PRF that can get a sweet spot
end

%% RX Simulated Samples
% build our pulse and map for pulses in slow time (N) and fast time (M) dimensions
function A = echoMatrix(p)
    t = (0:p.N-1).'/p.fs;
    tm = (0:p.M-1)*p.PRI;
    tau = 2*(p.R0 - p.v*tm)/p.c;
    A = p.A * lfm(t - tau, p) .* exp(-1j*2*pi*p.fc*tau); % received signal idealized with propagation roundtrip
end

%% TX LFM Build
% creates lfm pulse complex valued signal with phase shift
function s = lfm(t, p)
    s = exp(1j*pi*p.k*(t - p.Tp/2).^2) .* (t >= 0 & t < p.Tp);
    % (t - p.Tp/2).^2 centers the chirp, so it sweeps from −B/2 to +B/2
    % (t >= 0 & t < p.Tp) is the rect window: a logical mask that zeroes everything outside the pulse.

end

%% RawMatrix
% show measrued matrix with phase change and fast time axis

function rawMatrix(A, p)
    t = (0:p.N-1)/p.fs*1e6;
    m = 0:p.M-1;
    n0 = round((2*p.R0/p.c + p.Tp/2)*p.fs) + 1;
    
    figure
    subplot(2,1,1)
    imagesc(m, t, abs(A)); axis xy
    xlabel('Pulse m'); ylabel('Fast time (\mus)'); title('|A[n,m]|')
    subplot(2,1,2)
    plot(m, unwrap(angle(A(n0,:))), 'o', m, angle(A(n0,1)) + 2*pi*p.fD*p.PRI*m, '-')
    xlabel('Pulse m'); ylabel('Phase (rad)'); legend('measured', '\phi_0 + 2\pi f_D m T', 'Location', 'northwest')
end