function y_neu = RungeKuttaODE (y_start, incr, speed)
%This function computes the next step of iteration for the ODE. It uses the Runge-Kutta method of 4th order

k1 = RightSideODE(y_start, speed);
k2 = RightSideODE(y_start + 0.5*incr*k1, speed);
k3 = RightSideODE(y_start + 0.5*incr*k2, speed);
k4 = RightSideODE(y_start + incr*k3, speed);

 y_neu = y_start + incr*(k1+2*(k2 + k3)+k4)/6;
end