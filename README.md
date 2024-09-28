# Penalized Principal Machine (Penalized PM)
A repository for penalized Principal Machine solved by the Group Coordinate Descent algorithm (Penalized PM)

1. Functions below belong to folder 'R'
- fn_pplssvm.R: it includes the function 'pplssvm' for the penalized principal least squares SVM (P2LSM)
- fn_ppals.R: it includes the function 'ppalssvm' for the penalized principal asymmetric least square regression (P2AR)
- fn_ppl2svm.R: it includes the function 'ppl2svm' for the penalized principal L2-hinge least squares SVM (P2L2M)
- fn_pplr.R: it includes the function 'pplr' for the penalized principal logistic regression with GCD optimizer (P2LR)
- fn_wpplr.R: it includes the function 'wpplr' for the weighted penalized principal logistic regression (P2WLR)
- fn_spsdr.R: it includes the function 'spsvmSDR', which is a unified function for the penalized PM methods. It enjoys the proposed methods by changing the method argument.
- fn_test.R: it includes the example for testing function 'spsvmSDR'.
- fn_tune_pplssvm.R: it includes the function 'tune_pplssvm', which returns the optimal lambda for P2LSM via cross-validation
- fn_tune_ppalssvm.R: it includes the function 'tune_ppalssvm', which returns the optimal lambda for P2AR via cross-validation
- fn_tune_ppl2svm.R: it includes the function 'tune_ppl2svm', which returns the optimal lambda for P2L2M via cross-validation
- fn_tune_pplr.R: it includes the functions 'tune_pplr' and 'tune_wpplr', which returns the optimal lambda for P2(W)LR via cross-validation
- fn_minor_pPSDR.R: it includes the auxiliary functions for pPSDR (i.e., soft-thresholding, firm-thresholding operator, etc.).
- fn_penalized_logit_dc.R, fn_sparse_SIR.R, fn_tune_sparse_SIR.R: PPLR solved by DC algorithm, sparse SIR, and its tuning code, respectively.
2. The 'data' folder contains the Boston housing and Diagnostic Wisconsin Breast Cancer datasets.
3. The 'simulation' folder contains the replication R codes for the simulations and real data analysis.
  - fn_simulation_continuous.R: it includes the simulation code for Table 1.
  - fn_simulation_binary.R: it includes the simulation code for Table 2.
  - fn_simulation_time_n.R: it includes the simulation code for Figure 1.
  - fn_realdata_boston.R: it includes the simulation code for Table 3 and Figure 2.
  - fn_realdata_WDBC.R: it includes the simulation code for Figure 3.
  
