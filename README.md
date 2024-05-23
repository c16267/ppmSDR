# pPSDR
A repository for penalized Principal Sufficient Dimension Reduction with Group Coordinate Descent algorithm (pPSDR)

Below functions are belong to folder 'R'
- fn_pplssvm.R: it includes the function 'pplssvm' for the penalized principal least squares SVM (PPLS-SVM)
- fn_ppalssvm.R: it includes the function 'ppalssvm' for the penalized principal asymmetric least squares SVM (PPALS-SVM)
- fn_pplr.R: it includes the function 'pplr' for the penalized principal logistic regression with GCD optimizer (PPLR)
- fn_wpplr.R: it includes the function 'wpplr' for the weighted penalized principal logistic regression (wPPLR)
- fn_tune_pplssvm.R: it includes the function 'tune_pplssvm', which returns the optimal lambda for PPLS-SVM via CV
- fn_tune_ppalssvm.R: it includes the function 'tune_ppalssvm', which returns the optimal lambda for PPALS-SVM via CV
- fn_minor_pPSDR.R: it includes the auxiliary functions for pPSDR (i.e., soft-thresholding, firm-thresholding operator)

Below data sets are belong to folder 'data'
- wdbc.csv : Wisconsin Diagnostic Breast Cancer data
- BostonHousing.csv : Boston Housing data