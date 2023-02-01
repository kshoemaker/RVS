# RVS

Code and case study data for *"Bayesian feature selection for radiomics using reliability metrics"*, submitted to Frontiers in Genetics

## Running the Model

Running the model in both the simulation setting and the case study setting starts with the function `<setting>Harness`. 

### Simulation 

The function `SimHarness.m`, found [here](/MatlabCode), needs the following input variables:

* n - sample size
* balance - proportion of sample in group 1
* p1 - the number of discriminating features
* p2 - the number of noise features
* mu1 - the vector of true group mean values, 1 value per group
* informed - informed N or not?, 1 for informed, 0 for neutral

The harness function takes these inputs and simulates data according to those specifications, using the function `simdata.m`, also found [here](/MatlabCode). 

The code then calls the function `SimModelRun`, using the below values as arguments:

* Values produced by `simdata` or defined by the user
    * p1, p2 -  number of informative and noise variables 
    * X - train data
    * Xf - test data 
    * Y - train group membership
    * Yf - test group membership
    * N - reliablity parameter
* Values set in the Harness function that can be changed by the user if necessary
    * n_iter - number of iterations
    * bi - size of burn in
    * a, b - hyperparameters on sigma0j
    * ak, bk - hyperparameters on sigmaj1 and sigmaj2
    * alpha_0 and alpha_1 - parameters for probit prior
    * c - hyperparameter on sigma1 and sigma2
    * feature_thresh - threshold of ppi value for inclusion as variable
 
Output of Model: 

* numVar - number of selected variables
* numError - number of misclassified observations
* MargGam - PPI of each variable
* PostProb - The posterior probabilties for membership in each class
* mu01f, mu02f - Posterior mean of the variable coefficients
* mu01_gam, mu02_gam - iterations of mu, the variable coefficient
* GammaBI - iterations of Gamma, the latent selection parameter
* VS_AUC - computed AUC of variable selection
* ClassAUC - computed AUC of classification
* tpr_class - TPR for classification
* fpr_class - FPR for classification
* misclas - overall misclassification percentage
 
These variables are all saved in a `.mat` file, and the chains for the variable coefficient parameter, mu, are saved in group specific csvs. 

Example Usage: 

```
for i in 1:100
  SimHarness(iteration = i, n = 100, balance = 0.5, p1 = 4, p2 = 100, mu1 = [1, -1], informed = 0)
end 
```

### Case Study

#### Data Processing

The data should have the following variables:

* `X` - train data
* `Xf` - test data 
* `Y` - train group membership
* `Yf` - test group membership
* `N` - reliablity parameter

The set of possible variables should be centered and normality is assumed. In our case, where many of the variable were significantly non-normal, we applied the suggested power transform using BoxCox to the data.

The reliability parameter in our case is set by the data provided by R. Ger on the variablility of different radiomic features from machine to machine. The standard deviation is translated and inverted into a parameter that varies from 0 to 1 and indicates the reliability of the respective radiomic features. 

#### Calling the Model 

For the case study, `CaseStudyHarness` needs only one input variable, letting the code know whether you are using an informative reliablity parameter or not: 1 for informed, 0 for neutral

The code imports the processed data using a simple Matlab import function, `importCaseStudy.m`. 

The code then calls the function `CaseStudyModelRun`, using the below values as arguments:

* Values imported by `importCaseStudy`
    * X - train data
    * Xf - test data 
    * Y - train group membership
    * Yf - test group membership
    * N - reliablity parameter
* Values set in the Harness function that can be changed by the user if necessary
    * n_iter - number of iterations
    * bi - size of burn in
    * a, b - hyperparameters on sigma0j
    * ak, bk - hyperparameters on sigmaj1 and sigmaj2
    * alpha_0 and alpha_1 - parameters for probit prior
    * c - hyperparameter on sigma1 and sigma2
    * feature_thresh - threshold of ppi value for inclusion as variable

Output of MCMC: 

* numVar - number of selected variables
* numError - number of misclassified observations
* MargGam - PPI of each variable
* PostProb - The posterior probabilties for membership in each class
* mu01f, mu02f - Posterior mean of the variable coefficients
* mu01_gam, mu02_gam - iterations of mu, the variable coefficient
* GammaBI - iterations of Gamma, the latent selection parameter
* ClassAUC - computed AUC of classification
* tpr_class - TPR for classification
* fpr_class - FPR for classification
* misclas - overall misclassification percentage

These output variables are all saved in a `.mat` file, and the chains for the variable coefficient parameter, mu, are saved in group specific csvs. 

Note that this is similar to the simulation setting, with the exception of the variable selection AUC, as it is no longer possible to compute. 

Example usage: 
```
% where the data has been saved in `../CaseStudyData/`
CaseStudyHarness(1)
```

## MCMC Code 

In both cases, the `ModelRun` function calls the below functions to run the MCMC:
 
* `mainprogRR.m` initializes the class specific means and calls the variable selection code.
* `bvsRR.m` is the primary variable selection code, calculating the pieces of the marginal log likelihood and then calling the metropolis hastings code `metroRR.m`, then storing the results. 
* `metroRR.m` is the Metropolis Hastings step, containing the code for the add/delete/swap of the latent variable and updating the means, alphas, and intercepts as needed
  * `mmlRRdelta.m` and `mmlRgamma.m` compute the marginal log likelihoods of the means and the latent variables, respectively. 
  
 ## Data for our Case Study 
  
* Case study data is available [here](/CaseStudyData).
  * Data processing is also done in [HN_Aerts_Data.R](MatlabCode/HN_Aerts_Data.R)
* Reliability information is available [here](/CaseStudyData/LungProjectData). 
* Simulation data is created by the files beginning with `sim` in [MatLabCode](MatlabCode/)
