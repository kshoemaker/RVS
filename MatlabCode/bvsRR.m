function[mu_1_mat, mu_2_mat, ...
    log_prob, GammaA]= ...
bvsRR(X, Y, gam_prior, n_iter, r1, mu_j1, ...
    mu_j2, pPar, alpha_0, alpha_1, a, b, ak, bk, Q, d_k, ...
    aj, bj, n, p, mu_1_mat, mu_2_mat, N, h1)


% Randomly choose the starting gamma's to be "on"
if isempty(gam_prior)
    randp = randperm(p); 
    gam_prior = zeros(1, p); 
    onenum = randp(1:r1);  %choose the first r1 of the random indexes
    gam_prior(onenum) = 1; % set those indices equal to 1
end


disp(' ')
disp('------- RUNNING METROPOLIS ...')
disp(' ')

% --- Initialize parameters
% Counters and intial storage assignment
nPar = 1000; 
GammaA{1} = num2str(find(gam_prior));   
Uni(1) = r1; 
log_prob = zeros(1, n_iter + 1); 
ADgam = []; 
Sgam = []; 
counter = 1; 
counterD = 1;

% n1 is size of class 1, Y1 is members of class 1
n1 = sum(Y == 0); 
n2 = n-n1; 
Y1 = find(Y==0); 
Y2 = find(Y==1);

% varSel is the number of variables currently selected
varSel = sum(gam_prior); 

% QS is Q for only the selected variables
QS = Q(logical(gam_prior), logical(gam_prior));

% the mu parameter for only the selected variables
mu_1_sel = mu_j1(logical(gam_prior)); 
mu_2_sel = mu_j2(logical(gam_prior));


% ---- Calculate the marginal log-likelihood
% initialize zero vectors of length p
Avet = zeros(p, 1); 
B1vet = zeros(p, 1); 
B2vet = zeros(p, 1); 
D1vet = zeros(p, 1); 
D2vet = zeros(p, 1);

% set the a, ak and aj "prime"
a_p = a + n/2; 
ak_p1 = ak + n1/2; 
ak_p2 = ak + n2/2; 
aj_p = aj + 1/2;

% loop through variables, calculate the different pieces
%     of the marginal log likelihood
for i = 1:p
    Xj = X(:, i); 
    b_p = b + (Xj'*Xj)/2; % b prime 
    Xj1 = X(Y1, i);  % X of class 1 members and the current variable
    Xj2 = X(Y2, i);  % X of class 1 members and the current variable
    bk_p1 = bk + ((Xj1-ones(n1, 1) * mu_j1(i))'*(Xj1-ones(n1, 1) * ...
        mu_j1(i)))/2;  % bk prime for class 1
    bk_p2 = bk + ((Xj2-ones(n2, 1) * mu_j2(i))'*(Xj2-ones(n2, 1) * ...
        mu_j2(i)))/2;  % bk prime for class 2
    bj_p1 = bj + (mu_j1(i)^2)/2;  % bj twiddle for class 1
    bj_p2 = bj + (mu_j2(i)^2)/2;  % bj twiddle for class 2
    
    % pieces of the marginal log-likelihood
    % Note on all: we compute all of them, but only use the correctly indexed elements later

    % p(x_j(gamma^c) | gamma) 
    Avet(i) = -n*log(sqrt(2*pi)) + a*log(b) - a_p*log(b_p) ...
                + gammaln(a_p) - gammaln(a);
    % P(x_jk(gamma) | mu_0k, gamma), k = 1             
    B1vet(i) = -n1*log(sqrt(2*pi)) + ak*log(bk) - ak_p1*log(bk_p1) ...
                 - gammaln(ak) + gammaln(ak_p1);
    % P(x_jk(gamma) | mu_0k, gamma), k = 2
    B2vet(i) = -n2*log(sqrt(2*pi)) + ak*log(bk) - ak_p2*log(bk_p2) ...
                -gammaln(ak) + gammaln(ak_p2);
    % P(mu_jk(gamma^c) | gamma), k = 1
    D1vet(i) = -log(sqrt(2*pi)) + aj*log(bj) - aj_p*log(bj_p1) ...
                - gammaln(aj) + gammaln(aj_p);
    % P(mu_jk(gamma^c) | gamma), k = 2
    D2vet(i) = -log(sqrt(2*pi)) + aj*log(bj) - aj_p*log(bj_p2) ...
                - gammaln(aj) + gammaln(aj_p);
end

% p(mu_k | gamma) **combining both K!!!**
C = -varSel*log(pi) - varSel*log(2*h1) + d_k*log(det(QS)) ...
     + 2*gammaln((d_k + 1)/2) - 2*gammaln(d_k/2)...
    - (d_k+1) * (log(det(QS + (mu_1_sel*mu_1_sel')./(2*h1))) + log(det(QS + (mu_2_sel*mu_2_sel')./(2*h1))));



% combining the pieces, using the correctly indexed elements of each
loglik = C + sum(Avet(logical(gam_prior==0))) + ...
    sum(B1vet(logical(gam_prior))) + sum(B2vet(logical(gam_prior))) ...
    + sum(D1vet(logical(gam_prior==0))) ...
    + sum(D2vet(logical(gam_prior==0)));


% Prior on $\gamma$ (and z)
% alpha_1 = 6   % just setting the value of alpha 1

% aren't doing this anymore, we fixed alpha_1
%normrnd(w, t^2); % dray an init alpha_1
%z = mvnrnd(alpha_0 + alpha_1*N, ones(1, p)); % don't think this needs to be
%here? 
%prigam = sum(log(normcdf(alpha_0 + alpha_1*N)));


prigam_on = log(normcdf(alpha_0 + alpha_1*N));
prigam_off = log(1-normcdf(alpha_0 + alpha_1*N));
prigam = sum(prigam_on(logical(gam_prior))) + sum(prigam_off(logical(gam_prior==0)));


%prigam = mu*varSel+eta*gam_prior*MRF*gam_prior';

tic;
Uni = zeros(1, n_iter); 
prior_prob = loglik + prigam; % prior_prob is the sum of the two pieces
log_prob(counter) = prior_prob; % set the log_prob
NaccB = 0; 
acceptgam = 0;  
acceptmu1 = 0; 
acceptmu2 = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% calling the MH algorithm %%%%%%%%
%%%%%%% for n_iter iterations (set to 1000)

while counter<n_iter+1
    
%-------------------- Metropolis (nPar iterations)    
[UniSt,  ProbSt,  GammaSt,  loglik,  prior_prob, ...
    prigam,  mu_1_matSt, mu_2_matSt,  storeADgam,  ...
    storeSgam,  acceptgam, acceptmu1,  acceptmu2,  ...
    Avet,  B1vet,  B2vet,  C, D1vet, D2vet, NaccB] = ...
metroRR(X, Y1, Y2, gam_prior, loglik, prior_prob, ...
    prigam, mu_j1, mu_j2,pPar, alpha_0, alpha_1, ak, ...
    bk, Q, d_k, aj, bj, n1, n2, p, N, acceptgam, acceptmu1, ...
    acceptmu2, nPar, Avet, B1vet, B2vet, C, D1vet, D2vet, NaccB, h1);

%%%% Storing and printing the results of the MH 

log_prob(counter + 1:counter+nPar) = ProbSt; 
Uni(counter + 1:counter+nPar) = UniSt; 
GammaA = cat(1, GammaA, GammaSt); 
ADgam = cat(2,  ADgam,  (counter - 1) + storeADgam);  
Sgam = cat(2,  Sgam,  (counter - 1) + storeSgam);
counter = counter + nPar; 
gam_prior = zeros(1, p);
gam_prior(str2num(GammaSt{nPar})) = ...
    ones(1, length(str2num(GammaSt{nPar}))); 
mu_1_mat(:, (counter-nPar):(counter-1)) = mu_1_matSt;
mu_2_mat(:, (counter-nPar):(counter-1)) = mu_2_matSt;
mu_j1 = mu_1_matSt(:, end); mu_j2 = mu_2_matSt(:, end);

disp(sprintf('It.: %3d,  T: %3.1fs,  acpgam: %1.3f, acpB1: %1.2f, acpB2: %1.2f, Uni: %2d',counter-1,toc,(acceptgam/(counter-1)),(acceptmu1/NaccB),(acceptmu2/NaccB),Uni(counter-1))); %tic;  
tic;
end

% ------ Calculate acceptance ratios
acceptgam = acceptgam/n_iter; 
acceptmu1 = acceptmu1/NaccB; acceptmu2 = acceptmu2/NaccB;

disp(' ')
disp('--- End.') 
