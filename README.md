# pPSDR
A repository for penalized Principal Sufficient Dimension Reduction with Group Coordinate Descent algorithm (pPSDR)

Below functions are belong to folder 'R'
- fn_pplssvm.R: it includes the function 'pplssvm' for the penalized principal least squares SVM (PPLS-SVM)
- fn_ppals.R: it includes the function 'ppalssvm' for the penalized principal asymmetric least square regression (PPALS)
- fn_pplssvm.R: it includes the function 'ppl2svm' for the penalized principal L2-hinge least squares SVM (PPL2SVM)
- fn_pplr.R: it includes the function 'pplr' for the penalized principal logistic regression with GCD optimizer (PPLR)
- fn_wpplr.R: it includes the function 'wpplr' for the weighted penalized principal logistic regression (wPPLR)
- fn_tune_pplssvm.R: it includes the function 'tune_pplssvm', which returns the optimal lambda for PPLS-SVM via cross-validation
- fn_tune_ppalssvm.R: it includes the function 'tune_ppalssvm', which returns the optimal lambda for PPALS via cross-validation
- fn_tune_ppl2svm.R: it includes the function 'tune_ppl2svm', which returns the optimal lambda for PPL2-SVM via cross-validation
- fn_minor_pPSDR.R: it includes the auxiliary functions for pPSDR (i.e., soft-thresholding, firm-thresholding operator, and etc.)
