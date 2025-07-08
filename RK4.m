function [xNext] = RK4(x, t, dt, dotH, EarthParams, ChibisParams)
    k1 = ode(x, t, dotH, EarthParams, ChibisParams);
    k2 = ode(x + dt/2*k1, t + dt/2, dotH, EarthParams, ChibisParams);
    k3 = ode(x + dt/2*k2, t + dt/2, dotH, EarthParams, ChibisParams);
    k4 = ode(x + dt*k3, t + dt, dotH, EarthParams, ChibisParams);

    xNext = x + dt/6*(k1 + 2*k2 + 2*k3 + k4);
end