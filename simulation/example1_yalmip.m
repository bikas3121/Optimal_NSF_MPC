% example_yalmip
yalmip('clear')
clear 
clc 
close all

%% Plant - Reconstruction filter
Fs = 1e6;
Fc = 1e5;           % cutoff frequency
Wn = Fc/(Fs/2);
% Wn = pi/32;
n = 4;              % filter order

[b,a] = butter(n, Wn,"low");
[Ah, Bh, Ch, Dh] = tf2ss(b,a);

%% Optimization probelm setup
% Define variables
Pf = sdpvar(n);
Pg = sdpvar(n);
Wf = sdpvar(1,n);
Wg = sdpvar(n,1);
L = sdpvar(n);
mu_e  = sdpvar(1);


% Constraints
mu_eta = 1.5^2;
MA = [Ah*Pf + Bh*Wf , Ah; L ,Pg*Ah];
MB = [Bh ; Wg];
MC = [Ch*Pf + Dh*Wf ,  Ch];
MP = [Pf, eye(n); eye(n), Pg];
MC_tilde = [Wf, zeros(1,n)];

C1 = [MP MA MB; MA' MP zeros(2*n,1); MB' zeros(1,2*n) eye(1)];
C2 = [mu_e MC Dh'; MC' MP zeros(2*n,1); Dh zeros(1,2*n) eye(1)];
C3 = [mu_eta MC_tilde; MC_tilde' MP];

BC1 = eye(2*n + 2*n + 1 );
BC2 = eye(1 + 2*n + 1 );
BC3 = eye(1 + 2*n);

% Constraints
% F = [Pf >= eye(n), Pg>= eye(n), C1 >= BC1, C2 >= BC2, C3>= BC3];
F = [Pf >= 0, Pg>= 0, C1 >= 0, C2 >= 0, C3>= 0];

ops = sdpsettings('solver','sdpt3');
ops.verbose = 1;
ops.showprogress = 1;
ops.debug = 1;
optimize(F, mu_e, ops)

%% Extract values
Pf = value(Pf);
Pg = value(Pg);
Wf = value(Wf);
Wg = value(Wg);
L = value(L);

%%  Transformation matrices
inv_Pg = pinv(Pg);
Sf = Pf - inv_Pg;
P1 = [Pf, Sf; Sf, Sf]; 
P = pinv(P1);
U = [Pf, eye(n); Sf, zeros(n) ];

%% Optimal NTF
[Af, Bf, Cf, Df] = ntf(Ah, Bh, Pf, Pg, Wf, Wg, L); % state space
[bf, af] = ss2tf(Af, Bf, Cf, Df);           % transfer function

%% Noise shaping filter
NSF = tf(af-bf,af);
[A_nsf, B_nsf, C_nsf, D_nsf] = tf2ss(bf-af,af);

%% Closed loop system dynamics
Acl = [Ah, Bh*Cf; zeros(n) , Af];
Bcl = [Bh; Bf];
Ccl = [Ch, Dh*Cf];
Dcl = Dh;

%% 
MA1 = U'*P*Acl*U;
MB1 = U'*P*Bcl;
MC1 = Ccl*U;
MP1 = U'*P*U;

%% %% Plots
N = Fs/2;

[hf,wf] = freqz(bf,af, N, 'half', Fs);
[hlp,wlp] = freqz(b,a, N, 'half', Fs);
[hlp1,wlp1] = freqz(af,bf, N, 'half', Fs);
%NTF derived directly from butterworth filter
[hf1,wf1] = freqz(a,b, N, 'half', Fs);



sl = length(hf)/2;
figure(Name="Butterworth 2nd order")
plot(wlp(1:sl)/2*pi*1e-3, 20*log10(abs(hlp(1:sl))));
hold on 
plot(wf(1:sl)/2*pi*1e-3, 20*log10(abs(hf(1:sl))));
hold on 
plot(wf1(1:sl)/2*pi*1e-3, 20*log10(abs(hf1(1:sl))));
legend("$H(z)$", " $F_{opt}(z)$",  "$F_{mpc}(z)$", 'Interpreter','latex')
xlabel('Frequency (kHz)')
ylabel("Magnitude (dB)")
grid minor



figure
plot(wf(1:sl)/2*pi*1e-3, 20*log10(abs(hf(1:sl))));
hold on 
plot(wlp1(1:sl)/2*pi*1e-3, 20*log10(abs(hlp1(1:sl))));


legend("$F_{opt}(z)$","$H_{MPC}(z)$", 'Interpreter','latex')
grid minor
xlabel('Frequency (kHz)')
ylabel("Magnitude (dB)")


