function F = pend_F(xk,uk,params)
% F = pend_F(xk,uk,params)
% Jacobian of the state transition function for
% linearizing the pendulum system.
    K = params.K;
	m = params.m;
	L = params.L;
	g = params.g;
	dt = params.dt;
    F = [                 1        dt   0;
         -g*dt/L*cos(xk(1))  1-K*dt/m   0;
                          0         0   1];
end