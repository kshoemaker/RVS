function [ ] = ...
 SimHarness(iteration, n, balance, p1, p2, mu1, informed)

%%% This function writes out 3 files, one mat files with summary metrics, 
%%% two csv files with the chains of mu's for each group


%%% seed is set to the seed for this iteration
seeds = (1:100)+102389;
this_seed = seeds(iteration);
rng(this_seed) 

%% simulate data using function simdata

%%%% Input Variables (default provided for example)

    % n is sample size
    % n = 100;
    % balance is proportion of sample in group 1
    % balance = 0.5;
    % p1 is the number of discriminating features
    % p1 = 4;
    % p2 is the number of noise features
    % p2 = 100;
    % mu1 is the vector of true group mean values, 1 value per group
    % mu1 = [1, -1];
    % informed - informed N or not?, 1 for informed, 0 for neutral
    % informed = 1;

%%%% Output Variables

    % X - train data
    % Xf - test data 
    % Y - test group membership
    % Yf - train group membership
    % N - reliablity parameter

 [X, Xf, Y, Yf, N] = simdata(n, balance, p1, p2, mu1, informed);


%% run model using SimModelRun

%%%% Input Variables (default provided as example)

    % p1, p2 -  number of informative and noise variables 
        % (used to calculate performance) 
    % X - train data
    % Xf - test data 
    % Y - train group membership
    % Yf - test group membership
    % N - reliablity parameter
    % n_iter - number of iterations
         n_iter = 100000;
    % bi - size of burn in
         bi = 20000;
    % a, b - hyperparameters on sigma0j
         a = 3;
         b = 0.1;
    % ak, bk - hyperparameters on sigmaj1 and sigmaj2
         ak = 3;
         bk = 0.1;
    % alpha_0 and alpha_1 - parameters for probit prior
         alpha_0 = -2.75;
         alpha_1 = 3;
    % c - hyperparameter on sigma1 and sigma2
         c = 0.5;
    % feature_thresh - threshold of ppi value for inclusion as variable
         feature_thresh = 0.5;
    
%%%% Output Variables

    % numVar - number of selected variables
    % numError - number of misclassified observations
    % MargGam - PPI of each variable
    % PostProb - The posterior probabilties for membership in each class
    % mu01f, mu02f - Posterior mean of the variable coefficients
    % mu01_gam, mu02_gam - iterations of mu, the variable coefficient
    % GammaBI - iterations of Gamma, the latent selection parameter
    % VS_AUC - computed AUC of variable selection
    % ClassAUC - computed AUC of classification
    % tpr_class - TPR for classification
    % fpr_class - FPR for classification
    % misclas - overall misclassification percentage

[numVar, numError, MargGam, PostProb, mu01f, mu02f, mu01_gam, mu02_gam, GammaBI, VS_AUC, ClassAUC, tpr_class, fpr_class, misclas] = ...
SimModelRun(p1, p2, X, Xf, Y, Yf, N, n_iter, bi, a, b, ak, bk, alpha_0, alpha_1, c, feature_thresh);

if informed == 1
   type = "informative";
elseif informed == 0
    type = "noninformative";
end 

% create filenames
file_name = join(["./SimulationOutputFiles/iter", num2str(iteration), type, ".mat"], "");
csv1_name =  join(["./SimulationOutputFiles/iter", num2str(iteration), type, "chain1.csv"], "");
csv2_name =  join(["./SimulationOutputFiles/iter", num2str(iteration), type, "chain2.csv"], "");

% write out files
csvwrite(csv1_name, mu01_gam);
csvwrite(csv2_name, mu02_gam);
save(file_name, 'numVar', 'numError', 'MargGam', 'PostProb', 'mu01f', ...
    'mu02f', "GammaBI", 'VS_AUC', 'ClassAUC', 'tpr_class', 'fpr_class', 'misclas');
  

    



