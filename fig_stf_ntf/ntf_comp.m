clc
% close all
% clear all

Fs = 1022976;
Ts = 1/Fs;
z = tf('z',Ts);

tol = 1e-5;


% switch 1
% case 1
% FC1  = 10e4
    Fc1 = 10e4;
    Wn1 = Fc1/(Fs/2);
    [b1,a1] = butter(2,Wn1);
    H_ast1 = tf(b1,a1, Ts);
    
    F_ast1 = (H_ast1)*1/z;
    H1 = minreal(1/(1-F_ast1), tol);
    F1 = minreal((H1- 1)/H1, tol);
   
    NTF1 = minreal(1-F1, tol);
    STF1 = minreal(H1*NTF1, tol);
 
% FC2  = 20e4
    Fc2 = 20e4;
    Wn2 = Fc2/(Fs/2);
    [b2,a2] = butter(2,Wn2);
     H_ast2 = tf(b2,a2, Ts);
    F_ast2 = (H_ast2)*1/z;
    H2 = minreal(1/(1-F_ast2), tol);
    F2 = minreal((H2- 1)/H2, tol);
  
    NTF2 = minreal(1-F2, tol);
    STF2 = minreal(H2*NTF2, tol);

 % FC3  = 10e4
    Fc3 = 30e4;
    Wn3 = Fc3/(Fs/2);
    [b3,a3] = butter(2,Wn3);
    H_ast3 = tf(b3,a3, Ts);
    
    F_ast3 = (H_ast3)*1/z;
    H3 = minreal(1/(1-F_ast3), tol);
    F3 = minreal((H3- 1)/H3, tol);
   
    NTF3 = minreal(1-F3, tol);
    STF3 = minreal(H3*NTF3, tol);
 
% FC4  = 20e4
    Fc4 = 40e4;
    Wn4 = Fc4/(Fs/2);
    [b4,a4] = butter(2,Wn4);
     H_ast4 = tf(b4,a4, Ts);
    F_ast4 = (H_ast4)*1/z;
    H4 = minreal(1/(1-F_ast4), tol);
    F4 = minreal((H4- 1)/H4, tol);
  
    NTF4 = minreal(1-F4, tol);
    STF4 = minreal(H4*NTF4, tol);  
    
    % FC4  = 20e4
    Fc5 = 50e4;
    Wn5 = Fc5/(Fs/2);
    [b5,a5] = butter(2,Wn5);
     H_ast5 = tf(b5,a5, Ts);
    F_ast5 = (H_ast5)*1/z;
    H5 = minreal(1/(1-F_ast5), tol);
    F5 = minreal((H5- 1)/H5, tol);
  
    NTF5 = minreal(1-F5, tol);
    STF5 = minreal(H4*NTF5, tol);  
% end

% 
% figure
% h11 = bodeplot(STF1,{1e2,1e7});
% p11 = getoptions(h11);
% p11.FreqUnits = 'Hz';
% setoptions(h11, p11, 'PhaseVisible','off')
% 
% hold on 
% h12 = bodeplot(NTF1,{1e2,1e7});
% p12 = getoptions(h12);
% p12.FreqUnits = 'Hz';
% setoptions(h12, p12)    
% 
% hold on 
% h21 = bodeplot(STF2,{1e2,1e7});
% p21 = getoptions(h21);
% p21.FreqUnits = 'Hz';
% setoptions(h21, p21, 'PhaseVisible','off')
% 
% hold on 
% h22 = bodeplot(NTF2,{1e2,1e7});
% p22 = getoptions(h22);
% p22.FreqUnits = 'Hz';
% setoptions(h22, p22)    
% 
% hold on 
% h31 = bodeplot(STF3,{1e2,1e7});
% p31 = getoptions(h31);
% p31.FreqUnits = 'Hz';
% setoptions(h31, p31, 'PhaseVisible','off')
% 
% hold on 
% h32 = bodeplot(NTF3,{1e2,1e7});
% p32 = getoptions(h32);
% p32.FreqUnits = 'Hz';
% setoptions(h32, p32)    
% 
% hold on 
% h41 = bodeplot(STF4,{1e2,1e7});
% p41 = getoptions(h41);
% p41.FreqUnits = 'Hz';
% setoptions(h41, p41, 'PhaseVisible','off')
% 
% hold on 
% h42 = bodeplot(NTF4,{1e2,1e7});
% p42 = getoptions(h42);
% p42.FreqUnits = 'Hz';
% setoptions(h42, p42) 
% 
% legend('STF (100k)', 'NTF (100k)', 'STF (200k)', 'NTF (200k)', 'STF (300k)', 'NTF (300k)','STF (400k)', 'NTF (400k)')
% grid minor




figure 
h12 = bodeplot(NTF1,{1e2,1e7});
p12 = getoptions(h12);
p12.FreqUnits = 'Hz';
setoptions(h12, p12,'PhaseVisible','off')    


hold on 
h22 = bodeplot(NTF2,{1e2,1e7});
p22 = getoptions(h22);
p22.FreqUnits = 'Hz';
setoptions(h22, p22 ,'PhaseVisible','off')    



hold on 
h32 = bodeplot(NTF3,{1e2,1e7});
p32 = getoptions(h32);
p32.FreqUnits = 'Hz';
setoptions(h32, p32 ,'PhaseVisible','off')    



hold on 
h42 = bodeplot(NTF4,{1e2,1e7});
p42 = getoptions(h42);
p42.FreqUnits = 'Hz';
setoptions(h42, p42 ,'PhaseVisible','off') 

hold on 
h52 = bodeplot(NTF4,{1e2,1e7});
p52 = getoptions(h52);
p52.FreqUnits = 'Hz';
setoptions(h52, p52 ,'PhaseVisible','off') 

legend( 'NTF (100k)', 'NTF (200k)', 'NTF (300k)', 'NTF (400k)', 'NTF (500k)')
grid minor


%% freq resp vs ntf


figure 

h121 = bodeplot(H_ast1,{1e2,1e7});
p121 = getoptions(h121);
p121.FreqUnits = 'Hz';
setoptions(h121, p121,'PhaseVisible','off', 'ylim', [-90, 10])    


hold on 
h122 = bodeplot(NTF1,{1e2,1e7});
p122 = getoptions(h122);
p122.FreqUnits = 'Hz';
setoptions(h122, p122,'PhaseVisible','off', 'ylim', [-90, 10])    
hold on 

h321 = bodeplot(H_ast3,{1e2,1e7});
p321 = getoptions(h321);
p321.FreqUnits = 'Hz';
setoptions(h321, p321,'PhaseVisible','off', 'ylim', [-90, 10])   

hold on 
h322 = bodeplot(NTF3,{1e2,1e7});
p322 = getoptions(h322);
p322.FreqUnits = 'Hz';
setoptions(h322, p322 ,'PhaseVisible','off')    

hold on 
h521 = bodeplot(H_ast5,{1e2,1e7});
p521 = getoptions(h521);
p521.FreqUnits = 'Hz';
setoptions(h521, p521,'PhaseVisible','off', 'ylim', [-90, 10])   

hold on 
h522 = bodeplot(NTF5,{1e2,1e7});
p522 = getoptions(h522);
p522.FreqUnits = 'Hz';
setoptions(h522, p522 ,'PhaseVisible','off') 

legend('LPF (100k)', 'NTF (100k)', 'LPF (300k)',  'NTF (300k)','LPF (500k)', 'NTF (500k)')
grid minor



