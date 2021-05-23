clear all;
balancing_robot_params;
additional_params;

% x0 in the integrators in the electromechanical subsystem

% test 1
% x0_q = [0; pi/36]; 
% x0_qdot = [0; 0];
% ref = 0;
% disturbance = 0;

% test 2
% x0_q = [0; 0]; 
% x0_qdot = [0; 0];
% ref = 0.1/wheel.r;
% disturbance = 0;

% test 3
x0_q = [0; 0]; 
x0_qdot = [0; 0];
ref = 0.1/wheel.r;
disturbance = 5/drv.duty2V;

%% linearized model
% state q = [gamma, theta, dot_gamma, dot_theta]
% gamma = wheels, theta = body
M11 = M11q;
M12 = 2*gbox.N*(1-gbox.N)*mot.rot.Iyy+(body.m*body.zb+2*mot.rot.m*mot.rot.zb)*wheel.r;
M21 = M12;
M22 = M22q;
M = [M11,M12;...
     M21,M22];

G11 = 0;
G12 = 0;
G21 = 0;
G22 = -(body.m*body.zb+2*mot.rot.m*mot.rot.zb)*g;
G = [G11,G12;...
     G21,G22];
 
Minv = inv(M);
A21 = -Minv*G;
A22 = -Minv*Fv1;
A = [zeros(2,size(A21,2)),eye(size(A22,2));...
    A21,A22];
B = 2*gbox.N*mot.Kt/mot.R*[zeros(size(A,1)-size(Minv,1),size(Minv,2));Minv]*[1;-1];

sysC = ss(A,B,zeros(1,size(A,2)),zeros(1,size(B,2)));

%% exact discretization
Ts = 0.01;
sysD = c2d(sysC, Ts, 'zoh');
[Phi,Gamma,~,~] = ssdata(sysD);
H = [1,0,0,0];

%% feedforward
bb = zeros(size(Phi,1),1);
bb = [bb;ones(size(H,1))];
Nxu = [Phi - eye(size(Phi,1)), Gamma; H, zeros(size(H,1))]\bb;
Nx = Nxu(1:end-1);
Nu = Nxu(end);

%% Discrete LQR
% Bryson's rule
gammaBar = pi/18;
thetaBar = pi/360;
uBar = 1;

q11 = 1/gammaBar^2;
q22 = 1/thetaBar^2;
q33 = 0;
q44 = 0;
Q = diag([q11,q22,q33,q44]);

r11 = 1/uBar^2;
rho = 500;
R = r11*rho; 

K = dlqr(Phi,Gamma,Q,R);