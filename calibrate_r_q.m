function [r,q] = calibrate_r_q(S0,T,DiscFact,Fwd)
% calibrate r,q parameters of a Black/LV model 
%
% S0.. spot
% T.. vector of expiries
% DiscFact.. vector of discount factors
% Fwd.. vector of forwards

dt = diff([0 T]);
disc_fact = [1 DiscFact];
fwd = [S0 Fwd];

r = -1./dt.*log(disc_fact(2:end)./disc_fact(1:end-1));
q = -1./dt.*log(fwd(2:end)./fwd(1:end-1).*disc_fact(2:end)./disc_fact(1:end-1));

end
 