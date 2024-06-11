clc
% close all
% clear all

Fs = 1e6;
Ts = 1/Fs;
z = tf('z',Ts);
Fc = 49e4;
Wn = Fc/(Fs/2);
% Wn = 0.5;
[b,a] = butter(3,Wn,'low');
tol = 1e-5;
switch 1
case 1
    W = tf(b,a, Ts);
    F = (W)*1/z;
    C = minreal(1/(1-F), tol);
    Fb = minreal((C- 1)/C, tol);
    NTF = minreal(1/(1+ C*Fb), tol);
    STF = minreal(C*NTF, tol);
    inv_STF = minreal(1/STF);
    figure
    h1 = bodeplot(STF,{1e2,1e7});
    p1 = getoptions(h1);
    p1.FreqUnits = 'Hz';
    setoptions(h1, p1)
    hold on 
    h2 = bodeplot(NTF,{1e2,1e7});
    p2 = getoptions(h2);
    p2.FreqUnits = 'Hz';
    setoptions(h2, p2)    
    legend('STF','NTF')
    grid minor

case 2
    W = 1;
    F = W*1/z;
    C = 1/(1-F);
    Fb = (C- 1)/C;
    NTF = 1/(1+ C*Fb);
    STF = C*NTF;
   figure
    h1 = bodeplot(STF);
    p1 = getoptions(h1);
    p1.FreqUnits = 'Hz';
    setoptions(h1, p1)
    hold on 
    h2 = bodeplot(NTF);
    p2 = getoptions(h2);
    p2.FreqUnits = 'Hz';
    setoptions(h2, p2)    
    legend('STF','NTF')
    grid minor

case 3
    Wps = tf([2.245, -0.664],[1, 0.91], Ts);
    F = Wps*1/z;
    C = 1/(1-F);
    Fb = (C- 1)/C;
    NTF = minreal(1/(1+ C*Fb));
    STF = minreal(C*NTF);
 figure
    h1 = bodeplot(STF);
    p1 = getoptions(h1);
    p1.FreqUnits = 'Hz';
    setoptions(h1, p1)
    hold on 
    h2 = bodeplot(NTF);
    p2 = getoptions(h2);
    p2.FreqUnits = 'Hz';
    setoptions(h2, p2)    
    legend('STF','NTF')
    grid minor
end
%%







STF = STF; 
NTF = NTF;
b_stf = [STF.Numerator{1}];
a_stf = [STF.Denominator{1}];
b_ntf = [NTF.Numerator{1}];
a_ntf = [NTF.Denominator{1}];
n = Fs/100;
% 
% c_num = [C.num{1}];
% c_den = [C.den{1}];
% F_num = [Fb.num{1}];
% F_den = [Fb.den{1}];
% n = Fs/2;
% [h_c,w_c] = freqz(c_num,c_den, n, 'half', Fs);
% [h_f,w_f] = freqz(F_num,F_den, n, 'half', Fs);
%%
% 
% [h_stf,w_stf] = freqz(b_stf,a_stf, n, 'half', Fs);
% [h_ntf,w_ntf] = freqz(b_ntf,a_ntf, n, 'half', Fs);
% % 
% figure
% semilogx(w_stf, h_stf);
% hold on 
% semilogx(w_ntf, h_ntf);
% sl = 20000;
% figure
% plot(w_stf(1:sl)/2*pi*1e-3, 20*log10(abs(h_stf(1:sl))));
% hold on 
% plot(w_ntf(1:sl)/2*pi*1e-3, 20*log10(abs(h_ntf(1:sl))));
% % % hold on 
% % plot(w_c/2*pi*1e-3, 20*log10(abs(h_c)));
% % hold on 
% % plot(w_f/2*pi*1e-3, 20*log10(abs(h_f)));
% ylim ([-250, 50]);
% legend('STF','NTF')
% xlabel('Frequency (kHz)')
% ylabel('Magnitude (dB)')
% grid minor

