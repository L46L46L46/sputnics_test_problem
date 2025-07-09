function[controlConst] = controlParams(ChibisParams)
    t0 = ChibisParams.HMax / ChibisParams.MMax;
    OmegaMax = ChibisParams.omegaMax * t0;
    theta = ChibisParams.J(1, 1) / ChibisParams.J(3, 3); % в нашей постановке задачи J3 минимален
    % расчёт констант управления по материалам диссера Ткачева
    KOmega = - 2 * OmegaMax / ChibisParams.deltaMax * theta^2 / (2 * theta - 1) + ...
        + 1/ChibisParams.deltaMax * sqrt((2 * OmegaMax * theta^2 / (2 * theta - 1))^2 + ...
        + 4 * ChibisParams.HMax * t0 / ChibisParams.J(3, 3) * ChibisParams.deltaMax * theta^2 / (2 * theta - 1));
    KA = (2 * theta - 1) / (8 * theta^2) * KOmega^2;
    % возвращение к размерным величинам
    controlConst = struct();
    controlConst.komega = ChibisParams.J(3, 3) / t0 * KOmega;
    controlConst.ka = ChibisParams.J(3, 3) / t0^2 * KA;
end