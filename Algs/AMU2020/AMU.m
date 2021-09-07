function [H,errs,grads,ts] = AMU(A,options)
% AMU_SNMF
% solve \|A - G*G'\|_F^2
[maxiter,timelimit,~] = getopts(options);% parameters

G_0 = options.init; % initial G_0
F_0 = norm(A - G_0*G_0','fro')^2;
%fs = zeros(maxiter,1);
normA = sqrt(sum(sum(A.*A)));
ts = zeros(maxiter,1);
errs = zeros(maxiter,1);
grads = zeros(maxiter,1);

tic;
iter = 0;
iter_re = 0;
epsilon = 1.e-16;
while iter < maxiter
    %fs(iter+1) = 1/4*sqrt(F_0);
    ts(iter+1) = toc;
    H = G_0;
    AH = A*H; HH = H'*H;
    errs(iter+1) = sqrt(normA^2 + sum(HH(:).*HH(:)) -  2*sum(AH(:).*H(:)))/normA;
    ProH = H - max(H - (H*H'-A)*H,0);
    grads(iter+1) = max(abs(ProH(:))); 
    tic;
    %step2 update stage
    gamma = 1-3/(5+iter-iter_re);
    if iter == iter_re
        Y = G_0;
    else
        Y = max((1+gamma)*G_0 - gamma*G_t,epsilon);
    end
    G_new = Y.*((A*Y)./((Y*Y')*Y)).^(1/3);
    F_new = norm(A - G_new*G_new','fro')^2;
    if F_new>F_0
        iter_re = iter+1;
        G_new = G_t;
        F_new = F_0;
    end
    G_t = G_0;
    G_0 = G_new;
    F_0 = F_new;
    if  sum(ts)>=timelimit
        break;
    end
    iter = iter+1;
    
end
ts = ts(1:iter);
errs = errs(1:iter);
grads = grads(1:iter);
for i = 1:iter-1
    ts(i+1) = ts(i)+ts(i+1);
end
H = G_0;

end

function [maxiter,timelimit,ac] = getopts(options)
maxiter = 10000000; ac = 1.e-15;
if isfield(options,'maxiter'),  maxiter = options.maxiter; end
if isfield(options,'ac'),       ac = options.ac;             end
if isfield(options,'timelimit'),  timelimit = options.timelimit;   end

end
