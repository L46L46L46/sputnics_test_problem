function [dotH, LambdaRef] = flywhellControl(omega, omegaPrev, Q, LambdaRefPrevios, r, H, t, ChibisParams, controlConst, EarthParams)
    A = mQuat2dcm(Q, true);
    S = [A(2, 3) - A(3, 2);
         A(3, 1) - A(1, 3);
         A(1, 2) - A(2, 1)];
    rSSC = A * r(4:6);
    [omegaRef, LambdaRef] = RefOmegaQuat(r, LambdaRefPrevios, EarthParams);
    omegaRel = omega - A * omegaRef;

    gravityMoment = 3 * EarthParams.muE / norm(rSSC)^5 * cross(rSSC, ChibisParams.J * rSSC);
    externalMoment = gravityMoment;
    
    dotOmegaRef = A * (omegaRef - omegaPrev)/(t(2) - t(1));

    dotH = externalMoment + controlConst.ka * S + controlConst.komega * omegaRel - cross(omega, ChibisParams.J * omega) ...
        - ChibisParams.J * dotOmegaRef - cross(omega, H); 
end