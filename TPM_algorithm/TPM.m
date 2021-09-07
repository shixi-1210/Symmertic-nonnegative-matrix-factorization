function [H, errs, grads, ts, k] = TPM(A, options)

options_I.lambda = options.lambda;
options_I.timelimit  = options.timelimit;
options_I.ac = options.ac;
options_I.init = options.init; [H,errs_I,grads_I,ts_I] = Phase_I(A,options_I);

options_II.timelimit  = options.timelimit; options_II.t_old = ts_I(end);
options_II.ac = 1.e-7;
options_II.init = max(H,0); [H,errs_II,grads_II,ts_II] = Phase_II_LS(A,options_II);
errs  = [errs_I;errs_II]; grads = [grads_I;grads_II];  ts = [ts_I; ts_I(end) + ts_II];
k = length(errs_I);
