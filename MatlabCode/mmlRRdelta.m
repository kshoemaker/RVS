function [B1vetnew, B2vetnew] = mmlRRdelta(X, Y1, Y2, n1, n2, mu_j1, mu_j2, ak, bk, choice)

% ---- calculate the marginal log-likelihood
ak_p1 = ak + n1/2; 
ak_p2 = ak + n2/2; 


i = choice; 
Xj1 = X(Y1, i); 
Xj2 = X(Y2, i);
bk_p1 = bk + ((Xj1-ones(n1, 1)*mu_j1(i))'*(Xj1-ones(n1, 1)*mu_j1(i)))/2;
bk_p2 = bk + ((Xj2-ones(n2, 1)*mu_j2(i))'*(Xj2-ones(n2, 1)*mu_j2(i)))/2;

B1vetnew = -n1*log(sqrt(2*pi)) + ak*log(bk) - ak_p1*log(bk_p1) - gammaln(ak) + gammaln(ak_p1);
B2vetnew = -n2*log(sqrt(2*pi)) + ak*log(bk) - ak_p2*log(bk_p2) - gammaln(ak) + gammaln(ak_p2);
 