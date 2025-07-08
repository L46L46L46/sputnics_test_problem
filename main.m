clc
clear all
close all

dt = 0.1;
t = 0:dt:1e4;
N = length(t);

EarthParams = struct();
EarthParams.muE = 3.986e5; %в километрах
EarthParams.RE = 6371.302; %средний радиус в километрах
EarthParams.omegaE = 7.292 * 10^(-5);
EarthParams.J2 = 1052.8e-9;
EarthParams.mu0 = 7.8 * 1e6;
EarthParams.delta = 3/2 * EarthParams.mu0 * EarthParams.J2 * EarthParams.RE^2; % поправка на вторую гармионику

ChibisParams = struct();
ChibisParams.J = [3.0, -0.1, 0.0;
                 -0.1, 2.3, 0.2;
                 0.0, 0.2, 1.9];
ChibisParams.invJ = inv(ChibisParams.J);
ChibisParams.HMax = 0.2;
ChibisParams.MMax = 0.005;
% параметры для стабилизации - взяты из статьи МЮ
ChibisParams.deltaMax = 0.1;
ChibisParams.omegaMax = 1e-3;

% данные орбиты - взяты из статьи МЮ
inc = 56.7 * pi / 180;
hOrbit = 550;
JRot = 5e-4;

% задаём основные векторы
omega = zeros(3, N);
Q = zeros(4, N);
r = zeros(6, N);
H = zeros(3, N);

% расчёт постоянных управления
controlConst = controlParams(ChibisParams);

% задаём начальные данные

for k=1:N-1
    if k==1
        omegaPrev = [0; 0; 0];
        LambdaRefPrevios = [1; 0; 0; 0];
        Q(:, 1) = rand(4, 1);
        Q(:, 1) = Q(:, 1) / norm(Q(:, 1));
        omega(:, 1) = 1e-3 * ones(3, 1);
        random = 2 * pi * rand(1, 1);
        phi = pi/6;
        pOrbit = EarthParams.RE + hOrbit;
        r(1:3, 1) = sqrt(EarthParams.muE/pOrbit) * [cos(inc)*sin(phi); -cos(inc)*cos(phi); sin(inc)];
        r(4:6, 1) = pOrbit * [cos(phi); sin(phi); 0];
    
    else
        [omegaPrev, ~] = RefOmegaQuat(r(:, k), LambdaRefPrevios, EarthParams);
    end
    [dotH, LambdaRefPrevios] = flywhellControl(omega(:, k), omegaPrev, Q(:, k), LambdaRefPrevios, r(:, k), H(:, k), t, ChibisParams, controlConst, EarthParams);
    [omega(:, k+1), Q(:, k+1), r(:, k+1), H(:, k+1)] = integrator(omega(:, k), Q(:, k), r(:, k), H(:, k), dotH, t(k), dt, EarthParams, ChibisParams);
end

getPlot(t, r, 'r');
getPlot(t, omega, 'omega');
getPlot(t, Q, 'Q');