type = 'guassian';

switch type
    case 'syn'
        rng default;
        n = 200; r = 50;
        nexams = 100;
        
        N200R50 = cell(nexams,1);
        
        for i_exam = 1:nexams
            B = rand(n,r);
            E = randn(n,n);
            E = E + E';
            N200R50{i_exam,1} = {B,E};
        end
        
        save N200R50 N200R50;
    case 'guassian'
        Datasets = {...
            'PIE';...
            };
        n_datasets = length(Datasets);
        for nmlz = 2
            for i_data = 1:n_datasets
                data = Datasets{i_data,1};
                disp(data);
                if nmlz == 1
                    eval(['load Data/' data])
                    X = fea';
                else
                    eval(['load Data/' data ])
                    gnd = gnd - min(gnd) + 1;
                    X = fea;
                    X = X./sqrt(sum(X.*X,2));
                    X = X';
                end
                labels = gnd-min(gnd)+1;
                n = length(labels);
                n_class = length(unique(labels));
                D = Euclidean_dist(X);
                D2 = D.^2;
                q = floor(log2(n)) + 1;
                D_min = mink(D2,max(8,q),2);
                sigma = sqrt(D_min(:,8));
                sigma_q = D_min(:,q);
                A = exp(-bsxfun(@rdivide,D2,sigma*sigma'+eps));
                A = bsxfun(@times,A,D2<=max(sigma_q,sigma_q'));
                A = A-diag(diag(A));
                n = size(A,1);
                dk = sqrt(sum(A,2));
                A = bsxfun(@ldivide,dk+eps,A);
                A = bsxfun(@rdivide,A,dk'+eps);
                A = (A+A')/2;
                A = sparse(A);
                if nmlz == 1
                    eval(['save Data/A' data ' A labels'])
                else
                    eval(['save Data/A' data 'N2 A labels'])
                end
            end
        end
    case 'inner_product'
        Datasets = {...
            'TDT2';...
            };
        cc;
        load('TDT2.mat');
        X = fea;
        n = size(fea,1);
        X_unit = X./sqrt(sum(X.*X,2));
        W = X_unit*X_unit';
        p = floor(log2(n/max(gnd)))+1;
        b = maxk(W,p,2);
        A = W.*(W>=b(:,end));
        A = (A+A')/2;
        D = diag(1./sqrt(sum(A,2)));
        A = D*A*(D);
        labels = gnd;
end
