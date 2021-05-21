balrob_params;
ua2tau = 2*gbox.N*mot.Kt/mot.R*[1;...
                                -1];
l = body.zb;
mb = body.m;
Ibyy = body.Iyy;
w = 2*wheel.yb;
r = wheel.r;
mw = wheel.m;
Iwyy = wheel.Iyy;
zbrot = mot.rot.zb;
mrot = mot.rot.m;
Irotyy = mot.rot.Iyy;
n = gbox.N;
bw  = wheel.B;
bm  = mot.B;
bg  = gbox.B;
b   = n^2*bm+bg;

M11q = 2*Iwyy + 2*Irotyy*n^2 + (mb + 2*(mrot+mw))*r^2;
%M12q = 2*(1-n)*n*Irotyy + r*(l*mb + 2*mrot*zbrot)*cos(th);
%M21q = M12;
M22q = Ibyy + 2*(1-n)^2*Irotyy + mb*l^2 + 2*mrot*zbrot^2;
C11q = 0;
C22q = 0;
C21q = 0;
%C12q = -r*(mb*l + 2*mrot*zbrot)*sin(th)*dot_th;
Fv11 = 2*(b + bw);
Fv12 = -2*b;
Fv21 = Fv12;
Fv22 = 2*b;
Fv = [Fv11,Fv12;...
      Fv21,Fv22];
Fv1 = Fv + 2*n^2*mot.Kt*mot.Ke/mot.R*[1,-1;-1,1];
g1q = 0;
%g2q = -g*(mb*l + 2*mrot*zbrot)*sin(th);
%M22 = body.Iyy + 2*(1-gbox.N)^2*mot.rot.Iyy+body.m*