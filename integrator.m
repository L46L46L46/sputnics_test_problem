function [omegaNext, QNext, rNext] = integrator(omega, Q, r, t)
   x = [omega; Q; r];
   xNext = RK4(x, dt, EarthParams, ChibisParams);
   omegaNext = xNext(1:3);
   QNext = xNext(4:7);
   rNext(8:13) = x_new;
end