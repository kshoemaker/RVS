function [Cnew, Bvetnew, Dvetnew] = mmlRRbeta(gam_prior, X, Y1, Y2, n1, n2, p, Q, mu_j1new, mu_j2new...
        , ak, bk, aj, bj, d_k, l, Bvetnew, Dvetnew, k, h1)


varSel = sum(gam_prior); 
QS = Q(logical(gam_prior),  logical(gam_prior)); 
mu_1_sel = mu_j1new(logical(gam_prior)); 
mu_2_sel = mu_j2new(logical(gam_prior));
aj_p = aj + 1/2; 

if k==1
    ak_p1 = ak + n1/2; 
    Xj1 = X(Y1, l); 
    bj_p1 = bj+(mu_j1new(l)^2)/2;  
    bk_p1 = bk + ((Xj1-ones(n1, 1)*mu_j1new(l))'*(Xj1-ones(n1, 1)*mu_j1new(l)))/2;            
    Bvetnew(l) = -n1*log(sqrt(2*pi)) + ak*log(bk)-ak_p1*log(bk_p1)-gammaln(ak)+gammaln(ak_p1);         
    Dvetnew(l) = -log(sqrt(2*pi)) + aj*log(bj) - aj_p*log(bj_p1) ...
                - gammaln(aj) + gammaln(aj_p);

elseif k==2
    ak_p2 = ak + n2/2; 
    Xj2 = X(Y2, l);
    bj_p2 = bj+(mu_j2new(l)^2)/2;
    bk_p2 = bk + ((Xj2-ones(n2, 1)*mu_j2new(l))'*(Xj2-ones(n2, 1)*mu_j2new(l)))/2;
    Bvetnew(l) = -n2*log(sqrt(2*pi))+ak*log(bk)-ak_p2*log(bk_p2)-gammaln(ak)+gammaln(ak_p2);
    Dvetnew(l) = -log(sqrt(2*pi)) + aj*log(bj) - aj_p*log(bj_p2) ...
                - gammaln(aj) + gammaln(aj_p);

end

Cnew = -varSel*log(pi) - varSel*log(2*h1) + (d_k)*log(det(QS)) + ...
    2*(gammaln((d_k + 1)/2) - 2*gammaln((d_k)/2))...
    - (d_k+1) * (log(det(QS+(mu_1_sel*mu_1_sel')./(2*h1))) + log(det(QS+(mu_2_sel*mu_2_sel')./(2*h1))));
