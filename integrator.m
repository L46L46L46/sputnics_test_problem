function [] = integrator(satellite, k, t, dt, EarthParams, ChibisParams)
   x = [satellite.omega(:, k); satellite.Q(:, k); satellite.r(:, k); satellite.H(:, k)];
   xNext = RK4(x, t, dt, satellite.dotH(:, k), EarthParams, ChibisParams);
   satellite.omega(:, k+1) = xNext(1:3);
   satellite.Q(:, k+1) = xNext(4:7) / norm(xNext(4:7));
   satellite.r(:, k+1) = xNext(8:13);
   satellite.H(:, k+1) = xNext(14:16); 
end