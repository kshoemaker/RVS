function [ mu_1_mat, mu_2_mat, log_prob, GammaA] =  ... 
mainprogRR(X, Y, gam_prior, n_iter, r1, mu_j1, mu_j2, pPar, ... 
	alpha_0, alpha_1, a, b, ak, bk, Q, d_k, aj, bj, N, h1)


[n, p] = size(X);    % Number of samples and number of ROIs, respectively

% initialize the class specific means 
mu_1_mat = zeros(p, n_iter); 
mu_2_mat = zeros(p, n_iter);

disp(' ')
disp('------- Calling the variable selection code now  --------')
disp(' ')


[mu_1_mat, mu_2_mat, log_prob, GammaA]= ... 
bvsRR(X, Y, gam_prior, n_iter, r1, mu_j1, mu_j2, ...
	pPar, alpha_0, alpha_1, a, b, ak, bk, Q, d_k, aj, bj, n, ...
	p, mu_1_mat, mu_2_mat, N, h1);
   