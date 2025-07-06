function [xNext] = RK4(x, dt, EarthParams, ChibisParams)
    k1 = ode(x, u, t, EarthParams, ChibisParams);
    k2 = ode(x + dt/2*k1, u, t + dt/2, EarthParams, ChibisParams);
    k3 = ode(x + dt/2*k2, u, t + dt/2, EarthParams, ChibisParams);
    k4 = ode(x + dt*k3, u, t + dt, EarthParams, ChibisParams);

    xNext = x + dt/6*(k1 + 2*k2 + 2*k3 + k4);
end