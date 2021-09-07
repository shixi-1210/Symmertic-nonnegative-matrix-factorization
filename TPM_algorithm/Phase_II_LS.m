function [H, errs,grads,ts] = Phase_II_LS(A,options)

% Phase_II computes the Symmetric nonnegative matrix factorization based on
% projected gradient direction.
% min f(W)_{W>=0} =  0.25*\|WW^T-A\|_F^2
timelimit = options.timelimit;
[maxiter,ac,gamma,mu] = getopts(options);
t_old = options.t_old;
alpha_old = 1;
normA = norm(A,'fro');

% Compute f and gradient at the initial X
tic;
H = options.init;
AH = A*H; HH = H'*H;
f  = 0.25*norm(A-H*H','fro')^2;
df = H*HH - AH;

errs = ones(maxiter,1);
grads = ones(maxiter,1);
ts = zeros(maxiter,1);
ts(1) = toc;


% PGD itaration
for iter = 1:maxiter
    errs(iter) = sqrt(4*f)/normA;
    T = H - max(H-df,0);
    grads(iter) = max(abs(T(:)));
    tic;
    alpha = max(2*alpha_old,gamma); malpha = alpha;
    fa = f;  tau = 100; inner_iter = 0; 
    while inner_iter<20
        inner_iter = inner_iter + 1;
        Ha = H - alpha*df; HaP = max(Ha,0); dH = HaP-H; 
        norm_dfa = sum(dH(:).*df(:));
        fa  = 0.25*norm(A-HaP*HaP','fro')^2;
        if fa <= f + mu*norm_dfa
            break
        else
            c = -0.5*alpha*norm_dfa/(fa-f-norm_dfa);
            alpha = min(max(alpha/tau,c),alpha/10);
            malpha = min(malpha,alpha);
        end
    end
    
    f = fa;
    AHaP = A*HaP; HaPHaP = HaP'*HaP;
    df = HaP*HaPHaP - AHaP;
    H = HaP;
    alpha_old = alpha;
    ts(iter + 1) = toc;
    
    
    % Check convergence of NCG
    if  t_old + sum(ts)>timelimit
        break
    end
end
iter = iter + 1;
errs(iter) = sqrt(4*f)/normA;
T = H - max(H-df,0);
grads(iter) = max(abs(T(:)));
errs(iter+1:end) = [];
grads(iter+1:end) = [];
ts(iter+1:end) = [];
tsums = zeros(iter,1);
for ii = 1:iter
    tsums(ii) = sum(ts(1:ii));
end
ts = tsums;
end

function [maxiter,ac,gamma,mu] = getopts(options)
maxiter = 5000000; ac = 1.e-8;  gamma = 1.e-3; mu = 0.1;

if isfield(options,'maxiter'),  maxiter = options.maxiter; end
if isfield(options,'ac'),       ac = options.ac;             end
if isfield(options,'gamma'),    mu = options.gamma;             end
if isfield(options,'mu'),       mu = options.mu;             end

end