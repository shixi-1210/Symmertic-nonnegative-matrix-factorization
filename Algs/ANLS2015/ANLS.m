function [H,err,ts] = ANLS(A, k, params, max_E)
normA = norm(A,'fro');
timelimit = params.timelimit;
[n,r] = size(params.Hinit);
for i = 1:max_E
    
    [Ht,~,~,ts] = symnmf_anls(A, k, params);
    params.Hinit = rand(n,r);
    if i == 1
        H0 = Ht;
        err0 = norm(A-H0*H0','fro')/normA;
        t0 = ts;
        if err0 < 1.e-10 || ts(end) >= timelimit
            break;
        end
    end    
    if i>=2
        err = norm(A-Ht*Ht','fro')/normA;
        if err < 1.e-10 || ts(end) >= timelimit-t0(end)
            t0 = [t0; t0(end) + ts];
            err0 = err;
            break;
        end
        H0 = (err>=err0)*H0 + (err0>err)*Ht;
        err0 = min(err,err0);
        %这里变黄 是因为要拼接不得已emmmmm
        t0 = [t0; t0(end) + ts];
    end
end
H = H0;
ts = t0;
err = err0;
end