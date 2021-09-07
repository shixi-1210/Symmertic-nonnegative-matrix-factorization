function [H,errs,grads,ts] = Phase_I(A,options)

% Phase_I minimize  min f(H) =  0.25*||A - HH^T||_F^2 + lambda/2*\|(H)_-\|_F^2
timelimit = options.timelimit;
[lambda,rho,sigma,maxiter,ac] = getopts(options);
eta = 0.5*sigma/(sigma-rho);
ac_dfnorm = ac(1);
ac_steplength = ac(2);
normA = norm(A,'fro');
normA2 = normA^2;

tic;
% Compute f and gradient at the initial H
H = options.init; AH = A*H; HH = H'*H; Hn = min(H,0);
f  = 0.25*(normA2 + sum(HH(:).^2) -  2*sum(AH(:).*H(:))) + 0.5*lambda*sum(Hn(:).^2);
df = H*HH - AH + lambda*Hn;

errs = ones(maxiter,1);
grads = zeros(maxiter,1);
ts = zeros(maxiter,1);

D = -df;
df_norm = norm(df,'fro'); df_norm2 = df_norm^2;
df_d = -df_norm2;
ts(1) = toc;

% NCG itaration
for iter = 1:maxiter
    Hp = max(H,0);
    AHp = A*Hp; HpHp = Hp'*Hp; dfp = Hp*HpHp - AHp; 
    errs(iter) = sqrt((sum(HpHp(:).*HpHp(:)) -  2*sum(AHp(:).*Hp(:)) + normA2)/normA2);
    T = Hp - max(Hp-dfp,0);
    grads(iter) = max(abs(T(:)));

    tic;
    % 1. Choose a suitable step length c for the linear search
    rho_dfd = rho*df_d; sigma_dfd = sigma*df_d;
    a = 0; fa = f; dfa_d = df_d;
    
    % 1.1) Choose [a,b] with b not satisfying Wolfe condition (1)
    b = eta;
    while 1
        Hb = H + b*D;
        AHb = A*Hb; HbHb = Hb'*Hb;  Hbn = min(Hb,0);
        fb  = 0.25*(normA2 + sum(HbHb(:).^2) - 2*sum(AHb(:).*Hb(:))) + 0.5*lambda*sum(Hbn(:).^2);
        if fb > f + b*rho_dfd
            break;
        else
            b = 2*b;
        end
    end
    
    % 1.2) Determine c via shrinking the interval [a b]
    while 1
        % Choose c
        delta = b-a;
        c = a - 0.5*delta^2*dfa_d/(fb-fa-delta*dfa_d);
        c = max(c,eta*a+(1-eta)*b);
        Hc = H + c*D;
        AHc = A*Hc; HcHc = Hc'*Hc; Hcn = min(Hc,0);
        fc  = 0.25*(normA2 + sum(HcHc(:).^2) - 2*sum(AHc(:).*Hc(:))) + 0.5*lambda*sum(Hcn(:).^2);
        
        % Shrink the interval [a b] to [a c] or [c b]
        if fc > f + c*rho_dfd      % c does not satisfy Wolfe condition (1)
            b = c; fb = fc;
        else                       % c satisfy Wolfe conditions (1)
            dfc = Hc*HcHc - AHc + lambda*Hcn;
            dfc_d = sum(dfc(:).*D(:)); % = <D,dfc>
            if dfc_d < sigma_dfd   % c does not satisfy Wolfe condition (2)
                a = c; fa = fc; dfa_d = dfc_d;
            else                   % c satisfy Wolfe conditions (1, 2)
                break
            end
        end
        
        % Check the convergece of interval updating
        if  (b-a)/b<ac_steplength
            Hc = H;
            fc = f;
            dfc = df;
            break
        end
    end
    y = dfc-df;
    df_norm2_ = df_norm2;
    
    H = Hc; f = fc; df = dfc;
    df_norm = norm(df,'fro'); df_norm2 = df_norm^2;
    
    % Check convergence of NCG
    if max(abs(df(:))) < ac_dfnorm || sum(ts)>timelimit
        ts(iter + 1) = toc;
        break
    else
        % Update the conjugated gradient d
        beta = sum(y(:).*df(:))/df_norm2_;
        mu = 1.e-3;
        while 1
            Dn = -df + beta*D;
            d_norm = sqrt(sum(Dn(:).*Dn(:)));
            df_d = sum(df(:).*Dn(:));
            if df_d <mu*df_norm*d_norm
                break
            else
                beta = beta/2;
            end
        end
        D = Dn;
        ts(iter + 1) = toc;
    end
end
iter = iter + 1;
Hp = max(H,0);
AHp = A*Hp; HpHp = Hp'*Hp; dfp = Hp*HpHp - AHp;
errs(iter) = sqrt((sum(HpHp(:).*HpHp(:)) -  2*sum(AHp(:).*Hp(:)) + normA2)/normA2);
T = Hp - max(Hp-dfp,0);
grads(iter) = max(abs(T(:)));
errs(iter+1:end) = [];
%
grads(iter+1:end) = [];
ts(iter+1:end) = [];
tsums = zeros(iter,1);
for ii = 1:iter
    tsums(ii) = sum(ts(1:ii));
end
ts = tsums;
end

function [lambda,rho,sigma,maxiter,ac] = getopts(options)
lambda = 0.01; rho = 0.1; sigma = 0.4;  maxiter = 6000;
ac = [1.e-5, 1.e-15];

if isfield(options,'lambda'),      lambda = options.lambda;           end
if isfield(options,'rho'),        rho = options.rho;           end
if isfield(options,'sigma'),      sigma = options.sigma;       end
if isfield(options,'maxiter_NCG'),  maxiter = options.maxiter_NCG; end
if isfield(options,'ac'),       ac = options.ac;             end

end
