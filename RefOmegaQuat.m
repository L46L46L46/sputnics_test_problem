function [omegaRef, LambdaRef] = RefOmegaQuat(r, LambdaRefPrevios, EarthParams)
    % r = [v, r]
    % LambdaRefPrevios - прерыдущий кватернирон ориентации -
    % для сохранения знака

    %% базисные векторы ОСК

    j3 = r(4:6) / norm(r(4:6));
    j2 = cross(r(4:6), r(1:3)) / norm(cross(r(4:6), r(1:3)));
    j1 = cross(j2, j3) / norm(cross(j2, j3));

    R = [j1'; j2'; j3']; %матрица перехода из ИСК в ОпСК

    %возмущающее ускорение - от второй гармоники
    a = EarthParams.delta / norm(r(4:6))^5 * (r(4:6) * (5 * r(6)^2 / norm(r(4:6))^2 - 1) - 2 * [0; 0; r(6)]);

    j3Dot = (r(1:3) - j3 * dot(r(1:3), j3)) / norm(r(4:6));
    j2Dot = (cross(r(4:6), a) - j2 * dot(cross(r(4:6), a), j2)) / norm(cross(r(4:6), r(1:3)));
    j1Dot = cross(j2Dot, j3) + cross(j2, j3Dot);

    omegaOrb = [
            (dot(j2Dot, j3) - dot(j3Dot, j2)) / 2;
            (dot(j3Dot, j1) - dot(j1Dot, j3)) / 2;
            (dot(j1Dot, j2) - dot(j2Dot, j1)) / 2]; % опорная угловая скорость задана в Ref СК

    LambdaRef = mDcm2quat(R);

    %% сохранение знака
    if dot(LambdaRef(2:4), LambdaRefPrevios(2:4)) < 0
           LambdaRef = - LambdaRef;
    end
    omegaRef = R' * omegaOrb; % перевести в ИСК
end