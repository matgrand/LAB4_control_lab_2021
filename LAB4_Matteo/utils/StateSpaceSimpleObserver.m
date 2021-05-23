balrob_params;
fc = 0.35; 
Tc = 1/(2*pi*fc);

z = tf('z',Ts);
s = (z-1)/Ts; % FE
HFE = 1/(Tc*s+1);

s = (z-1)/Ts/z; % BE
HBE = 1/(Tc*s+1);
C = Ts/(Ts + Tc);

Hz = HFE;
Hcz = 1 - Hz;
N = 3;

Hw = (1-z^-N)/N/Ts;