function eqs = craneDAEFunction(t,in2,in3,param131,param132,param133,param134,param135,param136,param137,param138,param139)
%craneDAEFunction
%    EQS = craneDAEFunction(T,IN2,IN3,PARAM131,PARAM132,PARAM133,PARAM134,PARAM135,PARAM136,PARAM137,PARAM138,PARAM139)

%    This function was generated by the Symbolic Math Toolbox version 9.0.
%    25-Jan-2023 14:40:14

Dthetat = in2(5,:);
Dxt = in2(4,:);
Nc = in2(3,:);
YP81 = in3(1,:);
YP82 = in3(2,:);
YP84 = in3(4,:);
YP85 = in3(5,:);
theta = in2(2,:);
t2 = sign(Dxt);
t3 = sign(Nc);
t4 = cos(theta);
t5 = sin(theta);
t6 = param135+param136;
t7 = Dthetat.^2;
t8 = 1.0./t6;
eqs = [Nc.*1.0e+3-param134.*t6+param132.*param136.*(YP85.*t5+t4.*t7);YP85-(-param134.*t5+t4.*(t8.*(param131+param132.*param136.*t7.*(t5+param137.*t2.*t3.*t4))-param134.*param137.*t2.*t3)+(param132.^2.*param133.*param138.^2.*param139.*t7.*1.550313834014991e+1)./param136)./(param132.*(param136.*t4.*t8.*(t4-param137.*t2.*t3)-4.0./3.0));YP84+t8.*(-param131+param132.*param136.*(YP85.*t4-t5.*t7)+Nc.*param137.*t2.*t3.*1.0e+3);Dxt-YP81;Dthetat-YP82];