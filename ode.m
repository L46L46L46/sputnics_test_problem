%Функция принимает на вход вектор x, состоящий из абсолютной угловой
%скорости, кватерниона ориентации, скорости и положения аппарата, заданных
%в ИСК, параметры планеты, а так же параметры аппарата. Функция возвращает
%производную вектора x

function [xDot] = ode(x, EarthParams, ChibisParams)
    omega = x(1:3);
    Q = x(4:7);
    v = x(8:10);
    r = x(11:13);
    rSSC = mQuat2dcm(Q, true);

%% Определение механических моментов
    externalMoment = 3 * EarthParams.muE / norm(r) * cross(rSSC, ChibisParams.J * rSSC);
    controlMoment = -dotH - cross(omega, H);
    moment = externalMoment + controlMoment;

%% Расчёт производных
    dotOmega = ChibisParams.J * (moment - cross(omega, ChibisParams.J * omega));
    dotQ = 0.5 * [
            - dot(Q(2:4), omega);
            Q(1) * omega + cross(Q(2:4), omega)
            ];
    dotV = -EarthParams.muE / norm(r)^3 * r + ...
        + EarthParams.delta * r / norm(r)^5 * (5 * r(3)^2 / norm(r)^2 - 1) - ...
        - 2 * EarthParams.delta * [0; 0; r(3)];
    dotR = v;
    xDot = [dotOmega; dotQ; dotV; dotR];
end