classdef Satellite < handle

    properties
        omega;
        omegaRel;
        omegaRef;
        Q;
        QRel;
        QRef;
        r;
        H;
        dotH;
        omegaRotor;
        S;
        EarthParams;
    end

    methods

        function this = Satellite(hOrbit, inc, N, EarthParams)
            % основные матрицы
            this.omega = zeros(3, N);
            this.omegaRel = zeros(3, N);
            this.omegaRef = zeros(3, N);
            this.Q = zeros(4, N);
            this.QRel = zeros(4, N);
            this.QRef = zeros(4, N);
            this.r = zeros(6, N);
            this.H = zeros(3, N);
            this.dotH = zeros(3, N);
            this.omegaRotor = zeros(3, N);
            % орбитальные элементы
            phi = pi/6;
            pOrbit = EarthParams.RE + hOrbit;
            this.r(1:3, 1) = sqrt(EarthParams.muE/pOrbit) * [cos(inc)*sin(phi); -cos(inc)*cos(phi); sin(inc)];
            this.r(4:6, 1) = pOrbit * [cos(phi); sin(phi); 0];
            % угловые элементы
            this.Q(:, 1) = rand(4, 1);
            this.Q(:, 1) = this.Q(:, 1) / norm(this.Q(:, 1));
            this.omega(:, 1) = 1e-3 * ones(3, 1);
            % задание постоянных для удобства
            this.EarthParams = EarthParams;
            this.S = zeros(3, N);
        end

        function initRefRel(this, k)
            if (k == 1)
                [omegaRef, this.QRef(:, 1)] = RefOmegaQuat(this.r(:, 1), eye(4, 1), this.EarthParams);
                this.omegaRef(:, 1) = mQuat2dcm(this.Q(:, 1), true) * omegaRef;
                this.omegaRel(:, 1) = this.omega(:, 1) - this.omegaRef(:, 1);
                this.QRel(:, 1) = MNQuatMultiply(mQuatConj(this.QRef(:, 1)), this.Q(:, 1));
            else
                [omegaRef, this.QRef(:, k)] = RefOmegaQuat(this.r(:, k), this.QRef(:, k-1), this.EarthParams);
                this.omegaRef(:, k) = mQuat2dcm(this.Q(:, k), true) * omegaRef;
                this.omegaRel(:, k) = this.omega(:, k) - this.omegaRef(:, k);
                this.QRel(:, k) = MNQuatMultiply(mQuatConj(this.QRef(:, k)), this.Q(:, k));
            end
        end
    end
end
