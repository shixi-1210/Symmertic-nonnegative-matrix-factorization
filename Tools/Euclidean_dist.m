function [D,DD]=Euclidean_dist(A,B)
% Each column is a point.
[m,N] = size(A);
if nargin==1
    W = zeros(N);
    mm = 30000; is = 1;
    for i=1:ceil(m/mm)
        ie = min(m,i*mm);
        Y = double(A(is:ie,:));
        W = W+Y'*Y;
        is = ie+1;
        clear Y
    end
    a = diag(W);
    DD = abs(bsxfun(@plus,a,a') - 2*W);
else
    [ell,~] = size(B);
    if ell~=m 
        error(' The dimensions of A and B are not equal.');
    end
    a = sum(A.^2,1);
    b = sum(B.^2,1);
    DD = abs(bsxfun(@plus,a',b) - 2*A'*B);
end
D = [];
if nargout==1
	D = sqrt(DD);
end