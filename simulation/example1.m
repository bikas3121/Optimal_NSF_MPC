% Example 1

% Plant
Fs = 1e6;
n = 4;
[b,a] = butter(n, pi/32,"low");
[Ah, Bh, Ch, Dh] = tf2ss(b,a);

cvx_begin 
    variable Pf(n,n) symmetric
    variable Pg(n,n) symmetric
    variable Wf(1,n) 
    variable Wg(n,1)
    variable Wn(1,1)
    variable L(n,n)

    
cvx_end