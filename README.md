# Penalized Principal Machine (Penalized PM)
A repository for penalized Principal Machine solved by the Group Coordinate Descent algorithm (Penalized PM)

Below functions are belong to folder 'R'
- fn_pplssvm.R: it includes the function 'pplssvm' for the penalized principal least squares SVM (P2LSM)
- fn_ppals.R: it includes the function 'ppalssvm' for the penalized principal asymmetric least square regression (P2AR)
- fn_ppl2svm.R: it includes the function 'ppl2svm' for the penalized principal L2-hinge least squares SVM (P2L2M)
- fn_pplr.R: it includes the function 'pplr' for the penalized principal logistic regression with GCD optimizer (P2LR)
- fn_wpplr.R: it includes the function 'wpplr' for the weighted penalized principal logistic regression (P2WLR)
- fn_tune_pplssvm.R: it includes the function 'tune_pplssvm', which returns the optimal lambda for P2LSM via cross-validation
- fn_tune_ppalssvm.R: it includes the function 'tune_ppalssvm', which returns the optimal lambda for P2AR via cross-validation
- fn_tune_ppl2svm.R: it includes the function 'tune_ppl2svm', which returns the optimal lambda for P2L2M via cross-validation
- fn_tune_pplr.R: it includes the functions 'tune_pplr' and 'tune_wpplr', which returns the optimal lambda for P2(W)LR via cross-validation
- fn_minor_pPSDR.R: it includes the auxiliary functions for pPSDR (i.e., soft-thresholding, firm-thresholding operator, and etc.)
