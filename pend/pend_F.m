function F = pend_F(xak,uk,params)
% F = pend_F(xak,uk,params)
% Jacobian of the state transition function for
% linearizing the augmented pendulum system.
    K = params.K;
	m = params.m;
	L = params.L;
	g = params.g;
	dt = params.dt;

    F = [                  1        dt            0;
         -g*dt/L*cos(xak(1))  1-K*dt/m   dt/(m*L^2);
                           0         0            1];

end