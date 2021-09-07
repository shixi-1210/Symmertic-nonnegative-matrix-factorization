% Synetic data experiment code
% n_rep in line 37 control the repeat times of different initial points, in our experiment, we set n_rep = 20;
rng(2021);
i_iob = 1;
Jobs = {
    1,  'Syn', 'Syn', [];
    };
task = Jobs{i_iob,2};
s_data = Jobs{i_iob,3};
switch s_data
    case 'Syn'
        Datasets = {'N200R50';...
            };
    case 'RW'
        Datasets = {...
            'ACOIL20';...
            'AORL';...
            'APIE';...
            'ATDT2';...
            };
        nmlzs = [1,1,2,1];
        timelimits = [20 10 50 120];
end
Algs = {
    1, 'TPM';...
    2, 'ANLS2015';...
    3, 'BCD2016';...
    4, 'IBCD2017';...
    5, 'ANLS2018';...
    6, 'AMU2020';...
    };
n_algs = size(Algs,1);
switch task
    case 'Syn'
        for fakeloop = 1
            load N200R50.mat;
            n_rep = 1; %repeat times, in our experiment, we set n_rep = 20;
            ob_point = 20;
            Times = zeros(n_rep,n_algs,ob_point+1);
            Errs  = zeros(n_rep,n_algs,ob_point+1);
            Grads = zeros(n_rep,n_algs,ob_point+1);
            for i_rep = 1:n_rep
                H0 = N200R50{i_rep,1}{1,1};
                A = abs(H0*H0');
                n = size(A,1);
                n_class = size(H0,2);
                timelimit = 5;
                H_0 = rand(n, n_class);
                beta = sqrt(sum(sum((A*H_0).*H_0)))/norm(H_0'*H_0,'fro');
                H_0 = beta*H_0;
                for i_alg = 1:n_algs
                    disp(Algs{i_alg,2});
                    switch Algs{i_alg,2}
                        case 'TPM'
                            options_TPM.lambda = 10;
                            options_TPM.timelimit = timelimit;
                            options_TPM.init = H_0;
                            options_TPM.ac = [1.e-3, 1.e-10];
                            [H, errs, grads, ts] = TPM(A, options_TPM);
                            
                        case 'ANLS2015'
                            options_ANLS.timelimit = timelimit;
                            options_ANLS.Hinit = H_0;
                            [H, errs, grads, ts] = symnmf_anls(A,n_class,options_ANLS);
                        case 'BCD2016'
                            [H, errs, grads, ts] = SNMF_BCD_gill(A,50000,H_0',timelimit);
                            
                        case 'IBCD2017'
                            [H, errs, grads, ts] = SNMF_cyclic_vBSUM(A, 50000, H_0, timelimit);
                            
                        case 'ANLS2018'
                            [H,errs,grads,ts] =  HALS_sparse(A,H_0,H_0,50000,10,timelimit);
                            
                        case 'AMU2020'
                            options_AMU.timelimit = timelimit;
                            options_AMU.init = H_0;
                            [H,errs,grads,ts] = AMU(A,options_AMU);
                    end
                    if i_alg == 1
                        ts = reshape(ts,length(ts),1);
                        errs = reshape(errs,length(errs),1);
                        grads = reshape(grads,length(grads),1);
                        len_gap = length(ts)/ob_point;
                        ts = ts(1:floor(len_gap):end);
                        ts = [ts(1);ts(end-ob_point + 1:end)];
                        errs = errs(1:floor(len_gap):end);
                        errs = [errs(1); errs(end-ob_point + 1:end)];
                        grads = grads(1:floor(len_gap):end);
                        grads = [grads(1);grads(end-ob_point + 1:end)];
                        Times(i_rep,i_alg,:) = ts;
                        Errs(i_rep,i_alg,:) = errs;
                        Grads(i_rep,i_alg,:) = grads;
                        
                    elseif i_alg >1
                        ts = reshape(ts,length(ts),1);
                        errs = reshape(errs,length(errs),1);
                        grads = reshape(grads,length(grads),1);
                        len_gap = length(ts)/ob_point;
                        ts = ts(1:floor(len_gap):end);
                        ts = [ts(1);ts(end-ob_point + 1:end)];
                        errs = errs(1:floor(len_gap):end);
                        errs = [errs(1); errs(end-ob_point + 1:end)];
                        grads = grads(1:floor(len_gap):end);
                        grads = [grads(1);grads(end-ob_point + 1:end)];
                        Times(i_rep,i_alg,:) = ts;
                        Errs(i_rep,i_alg,:) = errs;
                        Grads(i_rep,i_alg,:) = grads;
                    end
                end
            end
            mt_tpm = mean(Times(:,1,:),1); mt_tpm = mt_tpm(:);
            me_tpm = mean(Errs(:,1,:),1); me_tpm = me_tpm(:);
            mg_tpm = mean(Grads(:,1,:),1); mg_tpm = mg_tpm(:);
            mt_anls = mean(Times(:,2,:),1); mt_anls = mt_anls(:);
            me_anls = mean(Errs(:,2,:),1); me_anls = me_anls(:);
            mg_anls = mean(Grads(:,2,:),1); mg_anls = mg_anls(:);
            mt_bcd = mean(Times(:,3,:),1); mt_bcd = mt_bcd(:);
            me_bcd = mean(Errs(:,3,:),1); me_bcd = me_bcd(:);
            mg_bcd = mean(Grads(:,3,:),1); mg_bcd = mg_bcd(:);
            mt_ibcd = mean(Times(:,4,:),1); mt_ibcd = mt_ibcd(:);
            me_ibcd = mean(Errs(:,4,:),1); me_ibcd = me_ibcd(:);
            mg_ibcd = mean(Grads(:,4,:),1); mg_ibcd = mg_ibcd(:);
            mt_hals = mean(Times(:,5,:),1); mt_hals = mt_hals(:);
            me_hals = mean(Errs(:,5,:),1); me_hals = me_hals(:);
            mg_hals = mean(Grads(:,5,:),1); mg_hals = mg_hals(:);
            mt_amu = mean(Times(:,6,:),1); mt_amu = mt_amu(:);
            me_amu = mean(Errs(:,6,:),1); me_amu = me_amu(:);
            mg_amu = mean(Grads(:,6,:),1); mg_amu = mg_amu(:);
            ft = 14;
            subplot(1,2,1);hold on;
            set(gca,'Fontsize',ft);
            set(gca,'Yscale','log');
            set(gca,'XGrid','on');
            set(gca,'YMinorTick','off');
            set(gca,'linewidth',1);
            set(gca,'position',[0.08 0.35 0.4 0.4*1.6]);
            plot(mt_tpm,me_tpm,'r-o','linewidth',1.5)
            plot(mt_anls,me_anls,'b->','linewidth',1.5)
            plot(mt_bcd,me_bcd,'c-*','linewidth',1.5)
            plot(mt_ibcd,me_ibcd,'k-p','linewidth',1.5)
            plot(mt_hals,me_hals,'m-d','linewidth',1.5)
            plot(mt_amu,me_amu,'g-<','linewidth',1.5)
            axis([0,timelimit,-inf,inf]);
            box on;
            legend('TPM','BCD','IBCD','ANLS','HALS','AMU','NumColumns',2,'fontsize',12);
            xlabel('CPU Time(s)','fontsize',ft);
            ylabel('E(t)','fontsize',ft);
            set(gca,'linewidth',1);
            subplot(1,2,2);hold on;
            set(gca,'Fontsize',ft);
            set(gca,'Yscale','log');
            set(gca,'XGrid','on');
            set(gca,'YMinorTick','off');
            set(gca,'linewidth',1);
            set(gca,'position',[0.55 0.35 0.4 0.4*1.6]);
            plot(mt_tpm,mg_tpm,'r-o','linewidth',1.5)
            plot(mt_anls,mg_anls,'b->','linewidth',1.5)
            plot(mt_bcd,mg_bcd,'c-*','linewidth',1.5)
            plot(mt_ibcd,mg_ibcd,'k-p','linewidth',1.5)
            plot(mt_hals,mg_hals,'m-d','linewidth',1.5)
            plot(mt_amu,mg_amu,'g-<','linewidth',1.5)
            axis([0,timelimit,-inf,inf]);
            box on;
            legend('TPM','BCD','IBCD','ANLS','HALS','AMU','NumColumns',2,'fontsize',12);
            xlabel('CPU Time(s)','fontsize',ft);
            ylabel('G(t)','fontsize',ft);
        end
end