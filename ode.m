%Функция принимает на вход вектор x, состоящий из абсолютной угловой
%скорости, кватерниона ориентации, скорости и положения аппарата, заданных
%в ИСК, параметры планеты, а так же параметры аппарата. Функция возвращает
%производную вектора x

function [xDot] = ode(x, t, controlExpected, EarthParams, ChibisParams)
    omega = x(1:3);
    Q = x(4:7);
    v = x(8:10);
    r = x(11:13);
    H = x(14:16);

    rSSC = mQuat2dcm(Q, true) * r;

%% Определение механических моментов
    gravityMoment = 3 * EarthParams.muE / norm(rSSC)^5 * cross(rSSC, ChibisParams.J * rSSC);
    externalMoment = gravityMoment;
    % реализуем управляющий момент в рамках маховиков
    controlMoment = zeros(3, 1);
    for i=1:3
        if (H(i) >= ChibisParams.HMax) % насыщение
            controlMoment(i) = 0;
        elseif (controlExpected(i) >= ChibisParams.MMax)
            controlMoment(i) = ChibisParams.MMax;
        else
            controlMoment(i) = controlExpected(i);
        end
    end
    moment = externalMoment + controlMoment;% + cross(omega, H);

%% Расчёт производных
    dotOmega = ChibisParams.invJ * (moment - cross(omega, ChibisParams.J * omega));
    dotQ = 0.5 * [
            - dot(Q(2:4), omega);
            Q(1) * omega + cross(Q(2:4), omega)
            ];
    dotV = -EarthParams.muE / norm(r)^3 * r + ...
        + EarthParams.delta * r / norm(r)^5 * (5 * r(3)^2 / norm(r)^2 - 1) ...
        - 2 * EarthParams.delta / norm(r)^5 * [0; 0; r(3)];
    dotR = v;
    dotH = - controlMoment - cross(omega, H);
    xDot = [dotOmega; dotQ; dotV; dotR; dotH];
end