function [controlExpected] = flywhellControl(satellite, t, k, Lambda, ChibisParams, controlConst, EarthParams)
    omega = satellite.omega(:, k);
    QStabilization = satellite.getQStabilization(Lambda, k);
    satellite.S(:, k) = 4 * QStabilization(1) * QStabilization(2:4);
    rSSC = mQuat2dcm(satellite.Q(:, k), true) * satellite.r(4:6, k);
    J = [ChibisParams.J(1, 1), 0, 0;
        0, ChibisParams.J(2, 2), 0;
        0, 0, ChibisParams.J(3, 3)]; % оставляем только диагональные элементы,
    % недиагональные дают возмущение, а диагональные знаем точно

    gravityMoment = 3 * EarthParams.muE / norm(rSSC)^5 * cross(rSSC, J * rSSC);
    externalMoment = gravityMoment;
    
    if k == 1
        dotOmegaRef = satellite.omegaRef(:, 1)/(t(2) - t(1));
    else
        dotOmegaRef = (satellite.omegaRef(:, k) - satellite.omegaRef(:, k-1))/(t(k) - t(k-1));
    end

    controlExpected = - externalMoment - controlConst.ka * satellite.S(:, k) - controlConst.komega * satellite.omegaRel(:, k) ...
        + cross(omega, J * omega) + J * dotOmegaRef;
end