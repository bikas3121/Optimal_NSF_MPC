% function to return NTF matrices:

function [Af, Bf, Cf, Df] = ntf(Ah, Bh, Pf, Pg, Wf, Wg, L)

inv_Pg = pinv(Pg);
Sf = Pf - inv_Pg;
inv_Sf = pinv(Sf);
Af = (Bh*Wf - inv_Pg*(L - Pg*Ah*Pf))*inv_Sf;
Bf = Bh - inv_Pg*Wg;
Cf = Wf*inv_Sf;
Df = 1;
end