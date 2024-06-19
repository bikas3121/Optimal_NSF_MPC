% Example 1

clear all
clc 
close all

%% 
Fs = 1e6;


%% Plant - Reconstruction filter
Fc = 1e4;           % cutoff frequency
% Wn = Fc/(Fs/2);
Wn = pi/32;
n = 4;              % filter order

[b,a] = butter(n, Wn,"low");
[Ah, Bh, Ch, Dh] = tf2ss(b,a);


%% Optimization 
cvx_clear
cvx_begin sdp
    variable Pf(n,n) symmetric
    variable Pg(n,n) symmetric
    variable Wf(1,n) 
    variable Wg(n,1)
    variable L(n,n)
    variable mu_e
    
    mu_eta = 1.5^2;
    
    MA = [Ah*Pf + Bh*Wf , Ah; L ,Pg*Ah];
    MB = [Bh ; Wg];
    MC = [Ch*Pf + Dh*Wf ,  Ch];
    MP = [Pf, eye(n); eye(n), Pg];
    MC_tilde = [Wf, zeros(1,n)];

    % const 1
    % C1 = [MP MA MB; MA' MP zeros(2*n,1); MB' zeros(1,2*n) eye(1)];
    % C2 = [mu_e MC Dh'; MC' MP zeros(2*n,1); Dh zeros(1,2*n) eye(1)];
    % C3 = [mu_eta MC_tilde; MC_tilde' MP];
    Pf >= eye(n);
    Pg >= eye(n);
    [MP MA MB; MA' MP zeros(2*n,1);
    MB' zeros(1,2*n) eye(1)] >= eye(2*n + 2*n + 1 );
    [mu_e MC Dh'; MC' MP zeros(2*n,1); Dh zeros(1,2*n) eye(1)] >= eye(1 + 2*n + 1 );
   [mu_eta MC_tilde; MC_tilde' MP] >= eye(1 + 2*n);
   minimize (mu_e)
cvx_end

%%  Transformation matrices
inv_Pg = pinv(Pg);
Sf = Pf - inv_Pg;
P1 = [Pf, Sf; Sf, Sf]; 
P = pinv(P1);
U = [Pf, eye(n); Sf, zeros(n) ];

%% Optimal NTF
[Af, Bf, Cf, Df] = ntf(Ah, Bh, Pf, Pg, Wf, Wg, L);

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

%%
[bf, af] = ss2tf(Af, Bf, Cf, Df);

% NTF is R(z) -1 thus
NSF = tf(bf-af,af);
[A_nsf, B_nsf, C_nsf, D_nsf] = tf2ss(bf-af,af);



%% Plots
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
legend("H(z)", " F_{opt}(z)",  "F_{mpc}(z)", intrepretor = 'latex')
grid minor
xlabel('Frequency (kHz)')
ylabel("Magnitude (dB)")

figure
plot(wf(1:sl)/2*pi*1e-3, 20*log10(abs(hf(1:sl))));
hold on 
plot(wlp1(1:sl)/2*pi*1e-3, 20*log10(abs(hlp1(1:sl))));


legend("F_{opt}(z)","H_{MPC}(z)", intrepretor = 'latex')
grid minor
xlabel('Frequency (kHz)')
ylabel("Magnitude (dB)")


