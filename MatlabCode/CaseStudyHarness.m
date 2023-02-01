function [ ] = ...
 CaseStudyHarness(informed)
%%% This function writes out 3 files, one mat file with summary metrics, 
%%% two csv files with the chains of mu's for each group

% informed - informed N or not?, 1 for informed, 0 for neutral
    % informed = 1;

%% Import the processed data

    % Data is from LungProjectData and CTPatientData
    % Processing is done in HN_Aerts_Data.R
    
importCaseStudy('../CaseStudyData/CaseStudyData_HPVOutcome.mat')


%% run model

%%%% Input Variables

    % X - train data
    % Xf - test data 
    % Y - test group membership
    % Yf - train group membership
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
        alpha_1 = 1;
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
    % ClassAUC - computed AUC of classification
    % tpr_class - TPR for classification
    % fpr_class - FPR for classification
    % misclas - overall misclassification percentage
    
[numVar, numError, MargGam, PostProb, mu01f, mu02f, mu01_gam, mu02_gam, ...
    GammaBI, ClassAUC, tpr_class, fpr_class, misclas] = ...
 CaseStudyModelRun(X, Xf, Y, Yf, N, n_iter, bi, a, b, ak, bk, alpha_0, ...
    alpha_1, c, feature_thresh);

if informed == 1
   type = "informative";
elseif informed == 0
    type = "noninformative";
end 

file_name = join(["./iter", num2str(iteration), type, ".mat"], "");
csv1_name =  join(["./iter", num2str(iteration), type, "chain1.csv"], "");
csv2_name =  join(["./iter", num2str(iteration), type, "chain2.csv"], "");

csvwrite(csv1_name, mu01_gam);
csvwrite(csv2_name, mu02_gam);
save(file_name, 'numVar', 'numError', 'MargGam', 'PostProb', 'mu01f', ...
    'mu02f', "GammaBI", 'VS_AUC', 'ClassAUC', 'tpr_class', 'fpr_class', ... 
        'misclas');
  

    



%% old, possibly useful, code

%real_saved_neut(s,:) = real;
%     poi_neut(s,:) = MargGam;
%     tpr_fpr_neut(s, :) = [tpr fpr] ;
%     misclas_neut(s) = numError; 
% 
%      [B fitInfo] = lasso(X,Y,'CV',10);
%      lambda1SE = fitInfo.Lambda1SE;
%      idxLambda1SE = fitInfo.Index1SE;
%  
%      coef = B(:,idxLambda1SE);
%      lasso_coef(s,:) = coef; 
%      coef0 = fitInfo.Intercept(idxLambda1SE);
%  
%      yhat = Xf*coef + coef0;
%      predYf = yhat > 0.5;
%      misclas_lasso(s) = 25- sum(predYf == Yf);
%         tpr = sum(predYf == 1 & Yf == 1) / sum(Yf == 1) ;
%         fpr = sum(predYf == 1 & Yf == 0) / sum(Yf == 0) ; 
%      tpr_fpr_lasso(s,:) = [tpr fpr]; 
     
     
% x = 1:1:104;
% y = MargGam;
% 
% X = [x; x];
% Y = [repelem(0,104); y.'];
% figure; hold on;
% plot(x,y,'*'),  xlim([1 104]),ylim([0 1.05])
% title('Var Selection');
% line(X,Y); 
% hold off;
 
% 
%  
% % [misclas misclas_neut misclas_lasso ]
% % 
%  tpr_vs = sum(sum(poi(:,1:4) > 0.5)) / (4*trials)
%  fpr_vs = sum(sum(poi(:,5:104) > 0.5)) / (100*trials)
% % 
%  tpr_neut_vs = sum(sum(poi_neut(:,1:4) > 0.5)) / (4*trials)
%  fpr_neut_vs = sum(sum(poi_neut(:,5:104) > 0.5)) / (100*trials)
% % 
%  tpr_lasso_vs = sum(sum(abs(lasso_coef(:,1:4)) > 1e-6)) / (4*trials)
%  fpr_lasso_vs = sum(sum(abs(lasso_coef(:,5:104)) > 1e-6)) /(100*trials)
% % 
% 
% [tpr_vs fpr_vs tpr_neut_vs fpr_neut_vs tpr_lasso_vs fpr_lasso_vs]
% 
%  [tpr_fpr (tpr_fpr(:,1) + (1-tpr_fpr(:,2)) - 1)]
%  
%  [tpr_fpr_neut (tpr_fpr_neut(:,1) + (1-tpr_fpr_neut(:,2)) - 1)]
%      
%  [tpr_fpr_lasso (tpr_fpr_lasso(:,1) + (1-tpr_fpr_lasso(:,2)) - 1)]
 

% %% matthew's correlation coeff
% 
%  tp = sum(sum(poi(:,1:4) > 0.5))
%  fp = sum(sum(poi(:,5:104) > 0.5))
% % tn = 100*trials - fp
% % fn = 4*trials - tp
% % 
% tpr = tp/(4*trials)
% fpr = fp/(100*trials)

% % mcc = (tp * tn - fp * fn) / (sqrt((tp+fp)*(tp + fn)*(tn+fp)*(tn + fn)))
% % 
%  tp = sum(sum(poi_neut(:,1:4) > 0.5))
%  fp = sum(sum(poi_neut(:,5:104) > 0.5))
% % tn = 100*trials - fp
% % fn = 4*trials - tp
% % 
% 
% tpr = tp/(4*trials)
% fpr = fp/(100*trials)
% % mcc_neut = (tp * tn - fp * fn) / (sqrt((tp+fp)*(tp + fn)*(tn+fp)*(tn + fn)))
% % 
%  tp = sum(sum(abs(lasso_coef(:,1:4)) > 1e-6))
%  fp = sum(sum(abs(lasso_coef(:,5:104)) > 1e-60))
% % tn = 100*trials - fp
% % fn = 4*trials - tp
% 
% tpr = tp/(4*trials)
% fpr = fp/(100*trials)
% % 
% % mcc_lasso = (tp * tn - fp * fn) / (sqrt((tp+fp)*(tp + fn)*(tn+fp)*(tn + fn)))
% % 
% % acc = (25*trials - sum(misclas)) / (25*trials)
% % acc_neut = (25*trials - sum(misclas_neut)) / (25*trials)
% % acc_lasso = (25*trials - sum(misclas_lasso)) / (25*trials)
% % % 
% % %% look at ROI of LASSO
% % 
% lambda = [10:-0.05:0.05, 0.04, 0.03, 0.02, 0.01, 0.005, 0.001, 0.0005, 0.0001, 0.00001, 0];
% [B fitInfo] = lasso(X,Y,'CV',10, 'Lambda', lambda);
% 
% n_lambda = size(lambda, 2);
% tpr_roc_lasso = zeros(n_lambda, 1);
% fpr_roc_lasso = zeros(n_lambda, 1);

% for beta_ind = 1:n_lambda 
%   var_sel = abs(B(:, beta_ind)) > 1e-6;
%   tpr_roc_lasso(beta_ind) = sum(var_sel(1:4)) / 4;
%   fpr_roc_lasso(beta_ind) = sum(var_sel(5:104)) / 100;
% end


%       
% % Fit lasso with 10-fold CV
% % Default setting for alpha is 1 (i.e. lasso)
% lasso_opts = glmnetSet;
% lambda = [10:-0.05:0.05, 0.04, 0.03, 0.02, 0.01, 0.005, 0.001, 0.0005, 0.0001, 0.00001, 0];
% lasso_cv = cvglmnet(X, Y,  'gaussian', 10, [], 'response', lasso_opts, 0);
%  
% % Results for all parameters tried
% beta_sel = lasso_cv.glmnet_object.beta;
% n_lambda = size(lambda, 2);
% tpr_roc_lasso = zeros(n_lambda, 1);
% fpr_roc_lasso = zeros(n_lambda, 1);
% for beta_ind = 1:n_lambda 
%   var_sel = beta_sel(:, beta_ind) ~= 0;
%   [tpr_roc_lasso(beta_ind), fpr_roc_lasso(beta_ind)] = tpr_fpr_var(var_sel, gamma_true);
% end
%  
% csvwrite(strcat('tpr_fpr_roc_lasso_model', num2str(model), '_iter', ...
%               num2str(cur_iter), '.csv'), [tpr_roc_lasso, fpr_roc_lasso]);


%  [tpr_fpr (tpr_fpr(:,1) + (1-tpr_fpr(:,2)) - 1)];
%  
%  [tpr_fpr_neut (tpr_fpr_neut(:,1) + (1-tpr_fpr_neut(:,2)) - 1)];
%      
%  [tpr_fpr_lasso (tpr_fpr_lasso(:,1) + (1-tpr_fpr_lasso(:,2)) - 1)];
%  
 
% clear
% seeds = (1:25)+10289;
% s = 3
% rng(seeds(s)*100)  
%  run('simdata_unbalanced.m')
%    % run('InputSimGBM.m')
%     run('InputCaseStudy.m')
% PostProb1 = PostProb(:,2)
%     
% rng(seeds(s)*100)  
% run('simdata_unbalanced_neutral.m')
% %run('InputSimGBM.m')
% run('InputCaseStudy.m')   
% PostProbNeut = PostProb(:,2)
% 
%  [B fitInfo] = lasso(X,Y,'CV',10);
%  lambda1SE = fitInfo.Lambda1SE;
%  idxLambda1SE = fitInfo.Index1SE;
% 
%  coef = B(:,idxLambda1SE);
%  lasso_coef(s,:) = coef; 
%  coef0 = fitInfo.Intercept(idxLambda1SE);
% 
%  yhat = Xf*coef + coef0;    
%  
% [PostProb1 PostProbNeut yhat Yf] 