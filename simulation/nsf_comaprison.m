%nsf comparision

Fs = 1e6;

%% NSF 1:
A1 = [0.0, 0.0; 1.0, 0.0];
B1 = [ 2.0; 0.0];
C1 = [1.0, -0.5];
D1 = [0.0];
[b_nsf1, a_nsf1] = ss2tf(A1,B1,C1,D1);
NSF1 = tf(b_nsf1, a_nsf1);
N = Fs/2;
[h_nsf1,w_nsf1] = freqz(b_nsf1,a_nsf1, N, 'half', Fs);


%% Optimal NSF
A2 =    [0.3531,   -0.3659,    0.0320; 0.9995,    0.0006,   -0.0008; 0.0007,    0.9987,    0.0013];
B2 = [ 1.0001; 0.0001; -0.0002];
C2 = [-1.4064    0.8164   -0.2450];
D2 = [1];
[b_nsf2, a_nsf2] = ss2tf(A2,B2,C2,D2);
NSF1 = tf(b_nsf2, a_nsf2);
N = Fs/2;
[h_nsf2,w_nsf2] = freqz(b_nsf2,a_nsf2, N, 'half', Fs);
%% Plots 
sl = length(h_nsf1);
figure()
plot(w_nsf1(1:sl)/2*pi*1e-3, 20*log10(abs(h_nsf1(1:sl))));
% hold on
% plot(w_nsf2(1:sl)/2*pi*1e-3, 20*log10(abs(h _nsf2(1:sl))));
% legend("NSF1","NSF Optimal")
xlabel("Frequency (kHz)")
ylabel("Magnitude (dB)")
grid minor