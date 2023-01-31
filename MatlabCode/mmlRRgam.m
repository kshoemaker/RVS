function [Cnew, ss] = mmlRRgam(gamnew, Q, mu_j1, mu_j2, d_k, alpha_0, alpha_1, N, h1)

varSel = sum(gamnew); 
QS = Q(logical(gamnew), logical(gamnew)); 
mu_1_sel = mu_j1(logical(gamnew)); 
mu_2_sel = mu_j2(logical(gamnew));

% calculating P(mu_jk(gamma) | gamma)
Cnew = -varSel*log(pi) - varSel*log(2*h1) + (d_k)*log(det(QS)) + ...
	2*(gammaln((d_k + 1)/2) - 2*gammaln((d_k)/2))...
    - (d_k+1) * (log(det(QS+(mu_1_sel*mu_1_sel')./(2*h1))) + log(det(QS+(mu_2_sel*mu_2_sel')./(2*h1))));

% Prior on $\gamma$
%ss = mu*varSel+eta*gamnew*MRF*gamnew';
%ss = sum(log(normcdf(alpha_0 + alpha_1*N)));

% update on 3.24.19 
% calculates P(gamma)
prigam_on = log(normcdf(alpha_0 + alpha_1*N));
prigam_off = log(1-normcdf(alpha_0 + alpha_1*N));
ss = sum(prigam_on(logical(gamnew))) + sum(prigam_off(~logical(gamnew)));

