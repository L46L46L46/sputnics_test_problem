function [omegaNext, QNext, rNext, HNext] = integrator(omega, Q, r, H, dotH, t, dt, EarthParams, ChibisParams)
   x = [omega; Q; r; H];
   xNext = RK4(x, t, dt, dotH, EarthParams, ChibisParams);
   omegaNext = xNext(1:3);
   QNext = xNext(4:7) / norm(xNext(4:7));
   rNext = xNext(8:13);
   HNext = xNext(14:16); 
end