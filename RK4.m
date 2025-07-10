function [xNext] = RK4(x, t, dt, controlExpected, EarthParams, ChibisParams)
    k1 = ode(x, t, controlExpected, EarthParams, ChibisParams);
    k2 = ode(x + dt/2*k1, t + dt/2, controlExpected, EarthParams, ChibisParams);
    k3 = ode(x + dt/2*k2, t + dt/2, controlExpected, EarthParams, ChibisParams);
    k4 = ode(x + dt*k3, t + dt, controlExpected, EarthParams, ChibisParams);

    xNext = x + dt/6*(k1 + 2*k2 + 2*k3 + k4);
end