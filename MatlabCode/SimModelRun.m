function [numVar, numError, MargGam, PostProb, mu01f, mu02f, mu01_gam, mu02_gam, GammaBI, VS_AUC, ClassAUC, tpr_class, fpr_class, misclas] = ...
SimModelRun(p1, p2, X, Xf, Y, Yf, N, n_iter, ...
         bi, a, b, ak, bk, ...
         alpha_0, alpha_1, c, feature_thresh)

% p1, p2 -  number of informative and noise variables 
    % (used to calculate performance) 
% X - train data
% Xf - test data 
% Y - train group membership
% Yf - test group membership
% N - reliablity parameter
% n_iter - number of iterations
    % n_iter = 100000;
% bi - size of burn in
    % bi = 20000;
% a, b - hyperparameters on sigma0j
    % a = 3;
    % b = 0.1;
% ak, bk - hyperparameters on sigmaj1 and sigmaj2
    % ak = 3;
    % bk = 0.1;
% alpha_0 and alpha_1 - parameters for probit prior
    % alpha_0 = -2.75;
    % alpha_1 = 3;
% c - hyperparameter on sigma1 and sigma2
    % c = 0.5;
% feature_thresh - threshold of ppi value for inclusion as variable
    % feature_thresh = 0.5;


%%  Setting up the Parameters and Hyperparameters %%

p = size(X, 2); % number of radiomic features

%MCMC setting
gam_prior = []; % selection latent variablel

% n_iter =  100000; % desired number of MCMC iterations
r1 = 2; % number of starting variables (gamma)
mu_j1 = zeros(p, 1); % class specific mu_jk
mu_j2 = zeros(p, 1); % class specific mu_jk

% Hyperparameters setting
% probability of Add/Delete vs Swap
pPar = 0.5; % as recommended in the motif code

% Prior on Sigma1 and Sigma2
% c = 0.6;   %  used for computing bj 
Q = c*eye(p); %% standard setting for the IW, times c
d_k = 3;      %% standard setting for the IW

% that implies for the Inv-Gamma Prior on sigmaj's
aj = d_k/2; 
bj = 2/c; 

% Prior on mu0's
h1 = 1; % scale parameter


%% Run 3 chains 
%%%%%% Calls the "main program" and runs the model %%%%%%

% Run the mainprog - Chain #1
[mu_1_mat, mu_2_mat, GammaA]= ... 
mainprogRR(X, Y, gam_prior, n_iter, r1, mu_j1, mu_j2, pPar, alpha_0, alpha_1, ...
    a, b, ak, bk, Q, d_k, aj, bj, N, h1);

% set the burn-in 
% bi = 20000;

% save the non-burnin sections
GammaBI1 = GammaA((bi+2):end);
mu_1_mat1 = mu_1_mat(:, (bi+1):end);
mu_2_mat1 = mu_2_mat(:, (bi+1):end);


clear GammaA   mu_1_mat mu_2_mat;
%save filename_1.mat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Run the mainprog - Chain #2
[mu_1_mat, mu_2_mat, GammaA]= ... 
mainprogRR(X, Y, gam_prior, n_iter, r1, mu_j1, mu_j2, pPar, alpha_0, alpha_1, a, b, ak, bk, Q, d_k, aj, bj, N, h1);


GammaBI2 = GammaA((bi+2):end);
mu_1_mat2 = mu_1_mat(:, (bi+1):end);
mu_2_mat2 = mu_2_mat(:, (bi+1):end);

clear GammaA  mu_1_mat mu_2_mat;
%save filename_2.mat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Run the mainprog - Chain #3
[mu_1_mat, mu_2_mat, GammaA]= ... 
mainprogRR(X, Y, gam_prior, n_iter, r1, mu_j1, mu_j2, pPar, alpha_0, alpha_1, a, b, ak, bk, Q, d_k, aj, bj, N, h1);


GammaBI3 = GammaA((bi+2):end);
mu_1_mat3 = mu_1_mat(:, (bi+1):end);
mu_2_mat3 = mu_2_mat(:, (bi+1):end);

clear GammaA  mu_1_mat mu_2_mat;
%save filename_3.mat


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Pool together 3 MCMC Chains

GammaA = [GammaBI1 ; GammaBI2; GammaBI3];
mu_1_mat = [mu_1_mat1 mu_1_mat2 mu_1_mat3];
mu_2_mat = [mu_2_mat1 mu_2_mat2 mu_2_mat3];


disp(' ')
disp('------- Calculating marginal probabilities for gamma (ROIs)')
disp(' ')
sss = size(GammaA); 
GammaBI=GammaA;  
aa = size(GammaBI);
div = 1; 
bb = round(aa(1)/div);
aaa = 1; 
freTot = zeros(p, 1);

for j=1:div
    bbb = bb*j;
    if j==div 
        bbb=size(GammaBI);
    end
    FreqUni = [];
    for i=aaa:bbb
        FreqUni = [FreqUni str2num(GammaBI{i})];
    end
    fre = tabulate(FreqUni); fsize = size(fre);
    m = p - fsize(1); Z0 = zeros(m, 1);
    fre = [fre(:,2);Z0]; freTot = fre+freTot;
    aaa = bbb+1;
end

disp(' ')
disp('------- Selected ROIs ...')
disp(' ')
MargGam = freTot./aa(1); 	
find(MargGam>feature_thresh)


    ROI_sel = find(MargGam>feature_thresh);
    mu01_PostMean = mean(mu_1_mat(ROI_sel, :), 2);
    mu02_PostMean = mean(mu_2_mat(ROI_sel, :), 2);

    test = isempty(ROI_sel);
    if test == 0
        mu01_gam = []; mu02_gam = [];
        GammaBI=GammaA; 
        aa = size(GammaBI);
        div = 1; 
        bb = round(aa(1)/div); 
        aaa = 1; 
        freTot = zeros(p, 1);
        bbb=size(GammaBI);
        mu01_gamBI = mu_1_mat(ROI_sel, :);
        mu02_gamBI = mu_2_mat(ROI_sel, :);
        for i=aaa:bbb
            if length(str2num(GammaBI{i}))==length(ROI_sel)
                if str2num(GammaBI{i})==ROI_sel'
                   
                    mu01_gam = [mu01_gam; mu01_gamBI(:, i)'];
                    mu02_gam = [mu02_gam; mu02_gamBI(:, i)'];
                end
            end    
        end
    end

    selected = size(ROI_sel);

 
%%    Prediction 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp(' ')
    disp('------- Prediction  -------')
    disp(' ')
    threshold = feature_thresh;
    numError = 0; 
    numVar = 0;

    biB = bi;
    Yf = double(Yf);
    Y1 = find(Y==0); 
    Y2 = find(Y==1); 
    Y1f = find(Yf==0); 
    Y2f = find(Yf==1);
    n1 = sum(Y==0); 
    n2 = sum(Y==1);
    n = n1 + n2;
    G=max(Y)+1;

    % a_k^prime
    ak_p1 = ak + n1/2; 
    ak_p2 = ak + n2/2;
    
    % a_k^new
    akf1 = ak_p1 + 1/2;
    akf2 = ak_p2 + 1/2;    

    numErrorV = [];
    numVarV = []; 


        gammaf = MargGam>threshold;
        gammaList = find(gammaf);
        Xgam = X(:, logical(gammaf));
        XgamC = X(:, logical(1-gammaf));
        Xfgam = Xf(:, logical(gammaf));
        XfgamC = Xf(:, logical(1-gammaf));
        Nk =[];
        X1gam = Xgam(Y1, :);
        X2gam = Xgam(Y2, :); 
        ROI_sel = gammaList;
        
        % mean of all iterations
         mu01_PostMean = mean(mu_1_mat(ROI_sel, :), 2);
         mu02_PostMean = mean(mu_2_mat(ROI_sel, :), 2);
        
        % mean of iterations of only when the selected variables are selected
         mu01f =  mean(mu01_gam); % Beta01_PostMean;
         mu02f =  mean(mu02_gam);
        
        % in cases where the selected group was never chosen, 
        % using the mean of each variable over all iterations
            % only happens when the above mean doesn't exist
        if isnan(mu01f) 
         mu01f =  mu01_PostMean;
         mu02f =  mu02_PostMean;
         mu01_gam = mu_1_mat(ROI_sel, :);
         mu02_gam = mu_2_mat(ROI_sel, :); 
        end 
         
        [nf,  ~] = size(Xf);
        ClassP = [];
        ClassPM = [];
        likelis = [];
        PpostM = [];
        PostProb = [];
        [nselected , ~] = size(gammaList);



     for i=1:nf        
            AA = 1; 
            BB = 1; 
            b1_p = 1; 
            b2_p = 1;
            b1f = 1; 
            b2f = 1;

            test = isempty(ROI_sel);
            if test == 0
                for j=1:nselected
                    Xj1 = X1gam(:, j);
                    Xj1f = Xfgam(i, j);

                    M = (Xj1-repmat(mu01f(j)', n1, 1));  
                    
                    %b_k^prime
                    b1_p = b1_p * (bk +  ((Xj1-repmat(mu01f(j)', n1, 1))'*(Xj1-repmat(mu01f(j)', n1, 1))/2));                                
                    %b_k^new
                    b1f = b1f * (bk + ((Xj1f-mu01f(j)')'*(Xj1f-mu01f(j)')+M'*M)/2);


                    Xj2 = X2gam(:, j);
                    Xj2f = Xfgam(i, j);

                    M = (Xj2-repmat(mu02f(j)', n2, 1)); 
                    
                    %b_k^prime
                    b2_p = b2_p * (bk + ((Xj2-repmat(mu02f(j)', n2, 1))'*(Xj2-repmat(mu02f(j)', n2, 1))/2));
                    %b_k^new
                    b2f = b2f * (bk + ((Xj2f-mu02f(j)')'*(Xj2f-mu02f(j)')+M'*M)/2);   
                    

                end
            end 

            AA =  exp(ak_p1 * log(b1_p) - akf1 * log(b1f)) ;  
            BB =  exp(ak_p2 * log(b2_p) - akf2 * log(b2f)) ;   

            likf1 = (-(1/2)*log((2*pi)) + log(n1/n) + log(AA) + gammaln(akf1) - gammaln(ak_p1));
            likf2 = (-(1/2)*log((2*pi)) + log(n2/n) + log(BB) + gammaln(akf2) - gammaln(ak_p2));
           
            
            Ppost = [likf1 likf2];
            likelis = [likelis ; Ppost] ;% storing the computed likelihoods of class 1 and 2
            ExpPpost = [exp(likf1) exp(likf2)];
            PpostM =  ExpPpost./sum(ExpPpost); 
            PostProb = [PostProb ; PpostM ]; % added to track the class probabilities
            ClassP = [ClassP find(Ppost==max(Ppost))]; % assigns class based on max likelihood
                                                       % is "correct" -
                                                       % matches Yf

     end    
    
    err = sum((ClassP-1-(Yf)')~=0);
    numError = err ;
    numVar = sum(gammaf);


%% Metrics 

% Classification TPR/FPR
misclas = numError/length(Yf);
tpr_class = sum((ClassP-1) == 1 & Yf' == 1) / sum(Yf' == 1) ;
fpr_class = sum((ClassP-1) == 1 & Yf' == 0) / sum(Yf' == 0) ;

% Classification AUC 
[~ , ~, ~, ClassAUC] = perfcurve(Yf, PostProb(:,2), 1);

% Variable Selection AUC
TrueVars = [zeros(1,p1) + 1 zeros(1,p2)];
[~, ~, ~, VS_AUC] = perfcurve(TrueVars, MargGam,1);


