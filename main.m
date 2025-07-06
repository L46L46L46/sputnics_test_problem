clc
clear all
close all

dt = 10;
t = 0:dt:1e5;
N = length(t);

EarthParams = struct();
EarthParams.muE = 3.986e5; %в километрах
EarthParams.RE = 6371.302; %средний радиус в километрах
EarthParams.omegaE = 7.292 * 10^(-5);
%EarthParams.mu0 = 7.8e6;
EarthParams.J2 = 1052.8e-9;
EarthParams.delta = 3/2 * EarthParams.mu0 * EarthParams.J2 * EarthParams.RE^2; % поправка на вторую гармионику

% данные орбиты -- взяты из статьи МЮ
inc = 56.7 * pi / 180;
hOrbit = 550;
HMax = 0.2;
MMax = 0.005;
JRot = 5e-4;

for k=1:N-1
    [omega(:, k+1), Q(:, k+1), r(k+1)] = integratior(omega(:, k), Q(:, k), r(:, k), t(k));
end