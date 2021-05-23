%% electromechanical
% input conversion
ua2tau1 = 2*gbox.N*mot.Kt/mot.R*[1; -1];

b = gbox.N^2*mot.B+gbox.B;

% inertia matrix
M11q = 2*wheel.Iyy + 2*mot.rot.Iyy*gbox.N^2 + (body.m + 2*(mot.rot.m+wheel.m))*wheel.r^2;
M22q = body.Iyy + 2*(1-gbox.N)^2*mot.rot.Iyy + body.m*body.zb^2 + 2*mot.rot.m*mot.rot.zb^2;

% non-linear functions of the Simulink model, u is the input:
% M12q -->  M21q = 2*(1-gbox.N)*gbox.N*mot.rot.Iyy + wheel.r*(body.zb*body.m + 2*mot.rot.m*mot.rot.zb)*cos(u);
% C12q -->  -wheel.r*(body.m*body.zb + 2*mot.rot.m*mot.rot.zb)*sin(u)*dot_th;
% g2q  -->  -g*(body.m*body.zb + 2*mot.rot.m*mot.rot.zb)*sin(u);

% viscous friction
Fv11 = 2*(b + wheel.B);
Fv12 = -2*b;
Fv21 = Fv12;
Fv22 = 2*b;
Fv = [Fv11, Fv12;...
      Fv21, Fv22];
Fv1 = Fv + 2*gbox.N^2*mot.Kt*mot.Ke/mot.R*[1,-1;-1,1];

% gravity
g1q = 0;

%% simple state -space observer
% discrete-time complementary filter to merge sensors
% cut-off frequency [Hz]
fc = 0.35;
Tc = 1/(2*pi*fc);

z = tf('z',Ts);

% Forward-Euler
s = (z-1)/Ts;
HFE = 1/(Tc*s+1);

% Backward-Euler
s = (z-1)/(Ts*z);
HBE = 1/(Tc*s+1);

% chosen low-pass
Hz = HFE;
Hcz = 1 - Hz;


% real-derivative filter to obtain angular speeds
N = 3;
Hw = (1-z^-N)/(N*Ts);

%% accelerometer
% non-linear functions of the Simulink model, u is the input:
% (formula 64)
% x-axis -->  wheel.r*u(5)*cos(u(2))+sens.mpu.zb*u(6)+g*sin(u(2))
% z-axis -->  wheel.r*u(5)*sin(u(2))-sens.mpu.zb*u(4)^2-g*cos(u(2))