clear;
Electromechanical;
StateSpaceSimpleObserver;
M11 = M11q;
M12 = 2*n*(1-n)*Irotyy+(mb*l+2*mrot*zbrot)*r;
M21 = M12;
M22 = M22q;
M = [M11,M12;...
     M21,M22];

G11 = 0;
G12 = 0;
G21 = 0;
G22 = -(mb*l+2*mrot*zbrot)*g;
G = [G11,G12;...
     G21,G22];
Minv = inv(M);
A21 = -Minv*G;
A22 = -Minv*Fv1;
A = [zeros(2,size(A21,2)),eye(size(A22,2));...
    A21,A22];
B = 2*n*mot.Kt/mot.R*[zeros(size(A,1)-size(Minv,1),size(Minv,2));Minv]*[1;-1];
sysC = ss(A,B,zeros(1,size(A,2)),zeros(1,size(B,2)));
sysD = c2d(sysC,Ts,'zoh');
[Phi,Gamma,~,~] = ssdata(sysD);
H = [1,0,0,0];

bb = zeros(size(Phi,1),1);
bb = [bb;ones(size(H,1))];
xx = [Phi - eye(size(Phi,1)),Gamma;H,zeros(size(H,1))]\bb;
Nx = xx(1:end-1);
Nu = xx(end);

gammaBar = pi/18;
omegaBar = pi/360;
uBar = 1;
q11 = 1/gammaBar^2;
q22 = 1/omegaBar^2;
q33 = 0;
q44 = 0;
Q = diag([q11,q22,q33,q44]);
r11 = 1/uBar^2;
rho = 500;
R = r11*rho; 
K = dlqr(Phi,Gamma,Q,R);