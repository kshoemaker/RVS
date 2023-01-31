
function [UniSt, ProbSt, GammaSt, loglik, prior_prob, prigam, mu_1_matSt, mu_2_matSt, storeADgam, ...
    storeSgam, acceptgam, acceptmu1, acceptmu2, Avet, B1vet, B2vet, C, D1vet, D2vet, NaccB]= ...
metroRR(X, Y1, Y2, gam_prior, loglik, prior_prob, prigam, mu_j1, mu_j2, pPar, alpha_0, alpha_1, ak, ...
bk, Q, d_k, aj, bj, n1, n2, p, N, acceptgam, acceptmu1, acceptmu2, nPar, Avet, B1vet, B2vet, C, D1vet, D2vet, NaccB, h1)

%%% we don't actually use these, I think. It stores the choices between
%%% add/delete and swap 
storeADgam=[]; 
storeSgam=[];

%%% initializing storage
GammaSt = cell(nPar, 1); 
ProbSt = zeros(1, nPar); 
UniSt = zeros(1, nPar); 
mu_1_matSt = zeros(p, nPar); 
mu_2_matSt = zeros(p, nPar);

for i=1:nPar
    %%%%%% M-H step for $\gamma$
    gamnew = gam_prior; 
    poszeros = find(gam_prior==0);

    % Forced to add variable if all are 0
    if length(poszeros)==p           
        chi='ADgam';
        choice=ceil(rand(1)*p);
        if choice==0
            choice=1;
        else
        end
        gamnew(choice)=1;        
        [Cnew, prigamnew] = mmlRRgam(gamnew, Q, mu_j1, mu_j2, d_k, alpha_0, alpha_1, N, h1);
        Adiff = -Avet(choice); 
        Ddiff = -D1vet(choice)-D2vet(choice);
        [B1vetnew, B2vetnew] = mmlRRdelta(X, Y1, Y2, n1, n2, mu_j1, mu_j2, ak, bk, choice);
        B1vet(choice) = B1vetnew; 
        B2vet(choice) = B2vetnew;
        Bdiff = B1vetnew+B2vetnew;
    elseif  isempty(poszeros)         % Forced to delete variable
        chi='ADgam';
        choice=ceil(rand(1)*p);
        if choice==0
            choice=1;
        else
        end
        gamnew(choice)=0;
        [Cnew, prigamnew] = mmlRRgam(gamnew, Q, mu_j1, mu_j2, d_k, alpha_0, alpha_1, N, h1);

        Adiff = Avet(choice);  
        Bdiff = -B1vet(choice)-B2vet(choice);
        Ddiff = D1vet(choice)+D2vet(choice);
    else            
        x=rand(1);
        if x<pPar
%-------------------Adding/Deleting
            chi='ADgam';
            choice=ceil(rand(1)*p);
            if choice==0
                choice=1;
            else
            end
            gamnew(choice) = abs(gam_prior(choice)-1);            
            [Cnew, prigamnew] = mmlRRgam(gamnew, Q, mu_j1, mu_j2, d_k, alpha_0, alpha_1, N, h1);
            % Adding
            if (gamnew(choice)-1)==0
                Adiff = -Avet(choice); 
                Ddiff = -D1vet(choice)-D2vet(choice);
                [B1vetnew, B2vetnew] = mmlRRdelta(X, Y1, Y2, n1, n2, mu_j1, mu_j2, ak, bk, choice);
                B1vet(choice) = B1vetnew; 
                B2vet(choice) = B2vetnew;
                Bdiff = B1vetnew+B2vetnew;
            % Deleting
            else
                Adiff = Avet(choice);  
                Bdiff = -B1vet(choice)-B2vet(choice);
                Ddiff = D1vet(choice)+D2vet(choice);
            end
        else
%-------------------Swapping
            chi='Sgam';
            poszeros = find(gam_prior==0);
            posones = find(gam_prior==1);
            choice1=ceil(rand(1)*length(poszeros));
            if choice1==0
                choice1 = 1;
            else
            end
            choice2 = ceil(rand(1)*length(posones));
            if choice1==0
                choice1 = 1;
            else
            end
            elem1=poszeros(choice1);    
            elem2=posones(choice2);     
            gamnew(elem1)=1;            
            gamnew(elem2)=0;            
            [Cnew, prigamnew] = mmlRRgam(gamnew, Q, mu_j1, mu_j2, d_k, alpha_0, alpha_1, N, h1);
            [B1vetnew, B2vetnew] = mmlRRdelta(X, Y1, Y2, n1, n2, mu_j1, mu_j2, ak, bk, elem1);
            B1vet(elem1) = B1vetnew; 
            B2vet(elem1) = B2vetnew;
            Adiff = Avet(elem2)-Avet(elem1); 
            Bdiff = B1vet(elem1)+B2vet(elem1)-B1vet(elem2)-B2vet(elem2);
            Ddiff = D1vet(elem2)+D2vet(elem2)-D1vet(elem1)-D2vet(elem1);
        end
    end
              
    move=min(0,  Adiff+Bdiff+Ddiff+Cnew-C+prigamnew-prigam );
    if move==0
        gam_prior = gamnew;     
        eval(['store', chi, '=[store', chi, ' i];'])        
        loglik = Cnew + sum(Avet(logical(gam_prior==0)))+sum(B1vet(logical(gam_prior)))+sum(B2vet(logical(gam_prior)))+sum(D1vet(logical(gam_prior==0)))+...
        sum(D2vet(logical(gam_prior==0)));
        prior_prob = loglik+prigamnew;                
        GammaSt{i} = num2str(find(gamnew));
        UniSt(i) = length(find(gamnew)); 
       	acceptgam = acceptgam+1;        
        prigam = prigamnew;
        C = Cnew;
        
    else if (rand(1) < exp(move))
            gam_prior = gamnew;     
            eval(['store', chi, '=[store', chi, ' i];'])        
            loglik = Cnew + sum(Avet(logical(gam_prior==0)))+sum(B1vet(logical(gam_prior)))+sum(B2vet(logical(gam_prior)))+sum(D1vet(logical(gam_prior==0)))+...
            sum(D2vet(logical(gam_prior==0)));
            prior_prob = loglik+prigamnew;
            GammaSt{i} = num2str(find(gamnew));
            UniSt(i) = length(find(gamnew)); 
            acceptgam = acceptgam+1;        
            prigam = prigamnew;
            C = Cnew;
        else
            GammaSt{i}=num2str(find(gam_prior));
            UniSt(i)=length(find(gam_prior));           
        end
    end


    
%%%%%%% Gibbs step for Alpha_1 and z
%%%%%%% we are fixing these for now. 

%------------ Update z_k

% for k=1:p
%     m = alpha_0 + alpha_1*N(k);
%     s = 1; 
%     if gamnew(k) == 1 
%         x_std = trandn((0-m)/s, inf);
%         z(k) = m + s*x_std;
%     elseif gamnew(k) == 0
%         x_std = trandn(-inf, (0-m)/s);
%         z(k) = m + s*x_std;
%     end
%     
% end
% 
% 
% %------------ Update alpha_1
% 
% mu_a = (sum(z - alpha_0*N)+w/t)/(sum(N.^2) + 1/t);
% v_a = (sum(N.^2) + 1/t)^(-1);
% alpha_1 = normrnd(mu_a, v_a);


%------------ Update intercepts $\mu_{01}$ and $\mu_{02}$
    varSel = sum(gam_prior); 
    varInd = find(gam_prior);
    sigBeta = 0.05;
    mu_j1new = mu_j1; 
    mu_j2new = mu_j2;
    for j=1:varSel
        NaccB = NaccB+1;
        l = varInd(j);
        % \mu_{01}
        k = 1;
        mu1_prior = mu_j1(l);
        mu1_new = mu1_prior + sigBeta.*randn(1);
        mu_j1new(l) = mu1_new;
        B1vetnew = B1vet; 
        D1vetnew = D1vet;
        [Cnew, B1vetnew, D1vetnew] = mmlRRbeta(gam_prior, X, Y1, Y2, n1, n2,p, Q, mu_j1new, mu_j2, ...
         ak, bk, aj, bj, d_k, l, B1vetnew, D1vetnew, k, h1); 
        logliknew = Cnew + sum(Avet(logical(gam_prior==0))) + sum(B1vetnew(logical(gam_prior))) + sum(B2vet(logical(gam_prior))) +...
        sum(D1vetnew(logical(gam_prior==0))) + sum(D2vet(logical(gam_prior==0)));        
        % if accepted, replace values with new ones 
        if (logliknew>loglik || logliknew-loglik > log(rand(1)))
            mu_j1(l) = mu1_new;
            acceptmu1 = acceptmu1 + 1;
            loglik = logliknew;
            prior_prob = logliknew + prigam;
            C = Cnew; 
            B1vet = B1vetnew; 
            D1vet = D1vetnew;
        end
        % \mu_{02}
        k = 2;
        mu2_prior = mu_j2(l);
        mu2_new = mu2_prior + sigBeta.*randn(1);
        mu_j2new(l) = mu2_new;
        B2vetnew = B2vet; 
        D2vetnew = D2vet;
        [Cnew, B2vetnew, D2vetnew] = mmlRRbeta(gam_prior, X, Y1, Y2, n1, n2, p, Q, mu_j1, mu_j2new...
        , ak, bk, aj, bj, d_k, l, B2vetnew, D2vetnew, k, h1);               
        logliknew = Cnew + sum(Avet(logical(gam_prior==0))) + sum(B1vet(logical(gam_prior))) + sum(B2vetnew(logical(gam_prior))) +...
        sum(D1vet(logical(gam_prior==0))) + sum(D2vetnew(logical(gam_prior==0)));        
        % if accepted, replace values with new ones
        if (logliknew>loglik || logliknew-loglik > log(randn(1))) 
                mu_j2(l) = mu2_new;
                acceptmu2 = acceptmu2 + 1;
                loglik = logliknew;
                prior_prob = logliknew+prigam;
                C = Cnew; B2vet = B2vetnew; D2vet = D2vetnew;
        end        
    end
ProbSt(i) = prior_prob;                
mu_1_matSt(:, i) = mu_j1;
mu_2_matSt(:, i) = mu_j2;
end
