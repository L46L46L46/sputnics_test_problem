clc
clear all
close all

dt = 10;
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
ChibisParams.JRotor = 5e-4;
% параметры для стабилизации - взяты из статьи МЮ
ChibisParams.deltaMax = 0.2;
ChibisParams.omegaMax = ChibisParams.HMax / ChibisParams.J(3, 3);

% данные орбиты - взяты из статьи МЮ
inc = 56.7 * pi / 180;
hOrbit = 550;

% задаём начальные данные и основные векторы
satellite = Satellite(hOrbit, inc, N, EarthParams);

% расчёт постоянных управления
controlConst = controlParams(ChibisParams);

for k=1:N-1
    satellite.initRefRel(k);
    satellite.controlExpected(:, k) = flywhellControl(satellite, t, k, ChibisParams, controlConst, EarthParams);
    satellite.omegaRotor(:, k) = satellite.H(:, k) / ChibisParams.JRotor;

    integrator(satellite, k, t(k), dt, EarthParams, ChibisParams);
    disp(100 * k / N);
end
satellite.initRefRel(N);

getPlot(t, satellite.omegaRel, 'omegaRel');
getPlot(t, satellite.QRel, 'QRel');
getPlot(t, satellite.omegaRotor, 'omegaRotor');
getPlot(t, satellite.S, 'S');
getPlot(t, satellite.H, 'H');
