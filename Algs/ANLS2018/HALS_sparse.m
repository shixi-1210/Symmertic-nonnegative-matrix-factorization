function [U,errs,grads,ts] = HALS_sparse(M,U,V,maxiter,lambda,timelimit,truelabel)

alpha = 0;
delta = 0;
if size(V,1) > size(V,2)
    V = V';
end

% Initialization
etime = tic; nM = norm(M,'fro')^2; 
[m,n] = size(M); [m,r] = size(U);
a = 0; errs = []; ts = []; grads = []; iter = 1; 

%initial error
errs = [errs norm(M-U*U','fro')/norm(M,'fro')];
ProH = U - max(U - (U*U'-M)*U,0);
grads = [grads max(abs(ProH(:)))]; 
ts = [ts 0];

if nargin >= 10
    [~,labelp] = max(V', [], 2);
    tempaccp = ClusteringMeasure(truelabel, labelp);
    acch(1) = tempaccp(1);
end


% % Scaling, p. 72 of the thesis
% eit1 = cputime; A = M*V'; B = V*V'; eit1 = cputime-eit1; j = 0;
% scaling = sum(sum(A.*U))/sum(sum( B.*(U'*U) )); U = U*scaling; 
% Main loop
while iter <= maxiter 
    % Update of U
    j = 1;
    if j == 1 % Do not recompute A and B at first pass
        % Use actual computational time instead of estimates rhoU
        M1 = [M sqrt(lambda)*V']; V1 = [V sqrt(lambda)*eye(r)];
        eit1 = cputime; A = M1*V1'; B = V1*V1'; eit1 = cputime-eit1; 
    end
    eit2 = cputime; eps = 1; eps0 = 1;
    U = HALSupdt(U',B',A',eit1,alpha,delta); U = U';
    % Update of V
    M2 = [M; sqrt(lambda)*U']; U2 = [U; sqrt(lambda)*eye(r)];
    eit1 = cputime; A = (U2'*M2); B = (U2'*U2); eit1 = cputime-eit1;
    eit2 = cputime; eps = 1; eps0 = 1; 
    V = HALSupdt(V,B,A,eit1,alpha,delta); 
    % Evaluation of the error e at time t
    if nargout >= 3
        errs = [errs norm(M-U*U','fro')/norm(M,'fro')];
        ProH = U - max(U - (U*U'-M)*U,0);
        grads = [grads max(abs(ProH(:)))]; 
        ts = [ts toc(etime)];
        if nargin >= 10
            [~,labelp] = max(V', [], 2);
            tempaccp = ClusteringMeasure(truelabel, labelp);
            acch(iter+2) = tempaccp(1);
        end
    end
    if grads(iter)<=1.e-7 || ts(iter)>=timelimit
        break
    end
    
    iter = iter + 1; j = 1; 
end
if size(V,1) < size(V,2)
    V = V';
end
ts = ts';
errs = errs';
grads = grads';




% Update of V <- HALS(M,U,V)
% i.e., optimizing min_{V >= 0} ||M-UV||_F^2 
% with an exact block-coordinate descent scheme
function V = HALSupdt(V,UtU,UtM,eit1,alpha,delta)
[r,n] = size(V); 
eit2 = cputime; % Use actual computational time instead of estimates rhoU
cnt = 1; % Enter the loop at least once
eps = 1; eps0 = 1; eit3 = 0;
while cnt == 1 || (cputime-eit2 < (eit1+eit3)*alpha && eps >= (delta)^2*eps0)
    nodelta = 0; if cnt == 1, eit3 = cputime; end
        for k = 1 : r
            deltaV = max((UtM(k,:)-UtU(k,:)*V)/UtU(k,k),-V(k,:));
            V(k,:) = V(k,:) + deltaV;
            nodelta = nodelta + deltaV*deltaV'; % used to compute norm(V0-V,'fro')^2;
            if V(k,:) == 0, V(k,:) = 1e-16*max(V(:)); end % safety procedure
        end
    if cnt == 1
        eps0 = nodelta; 
        eit3 = cputime-eit3; 
    end
    eps = nodelta; cnt = 0; 
end