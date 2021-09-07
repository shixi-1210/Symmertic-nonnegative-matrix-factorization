% Synetic experiment
rng(2021);
i_iob = 1;
Jobs = {
    1,  'RW' , 'RW', []
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
n_datasets = length(Datasets);
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
    case 'RW'
        for fakeloop = 1
            n_rep = 1;
            ob_point = 10;
            Clus = zeros(n_datasets,n_rep,n_algs,2);
            for i_data = 1:n_datasets
                nmlz = nmlzs(i_data);
                data = Datasets{i_data,1};
                disp(data);
                if nmlz == 1
                    eval(['load Data/' data])
                    gnd = labels;
                else
                    eval(['load Data/' data 'N'])
                    gnd = labels - min(labels) + 1;
                end
                n = length(labels);
                n_class = length(unique(labels));
                
                timelimit = timelimits(i_data);
                Times = zeros(n_rep,n_algs,ob_point+1);
                Errs  = zeros(n_rep,n_algs,ob_point+1);
                Grads = zeros(n_rep,n_algs,ob_point+1);
                for i_rep = 1:n_rep
                    H_0 = rand(n, n_class);
                    beta = sqrt(sum(sum((A*H_0).*H_0)))/norm(H_0'*H_0,'fro');
                    H_0 = beta*H_0;
                    for i_alg =1:n_algs
                        disp(Algs{i_alg,2});
                        switch Algs{i_alg,2}
                            case 'TPM'
                                options_TPM.lambda = .01;
                                options_TPM.timelimit = timelimit;
                                options_TPM.init =  H_0;
                                options_TPM.ac = [1.e-4, 1.e-10];
                                [H, errs, grads, ts] = TPM(A, options_TPM);
                                
                            case 'ANLS2015'
                                options_ANLS.timelimit = timelimit;
                                options_ANLS.Hinit = H_0;
                                [H, errs, grads, ts] = symnmf_anls(A,n_class,options_ANLS);
                                
                            case 'BCD2016'
                                [H, errs, grads, ts] = SNMF_BCD_gill(A,50000,H_0',timelimit);
                                H = H';
                            case 'IBCD2017'
                                [H, errs, grads, ts] = SNMF_cyclic_vBSUM(A, 50000, H_0, timelimit);
                                
                            case 'ANLS2018'
                                [H,errs,grads,ts] =  HALS_sparse(A,H_0,H_0,50000,10,timelimit);
                                
                            case 'AMU2020'
                                options_AMU.timelimit = timelimit;
                                options_AMU.init = H_0;
                                [H,errs,grads,ts] = AMU(A,options_AMU);
                        end
                        [~,idx] = max(H,[],2);
                        Clus(i_data,i_rep,i_alg,1) = calAC(idx,gnd);
                        Clus(i_data,i_rep,i_alg,2) = calMI(idx,gnd);
                        
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
                subplot(2,2,i_data);hold on;
                if i_data == 1
                    subplot('Position',[0.15 0.59 0.25 0.25*1.4]);
                elseif i_data == 2
                    subplot('Position',[0.47 0.59 0.25 0.25*1.4]);
                elseif i_data == 3
                    subplot('Position',[0.15 0.09  0.25 0.25*1.4]);
                elseif i_data == 4
                    subplot('Position',[0.47 0.09  0.25 0.25*1.4]);
                end
                hold on;
                ft = 13;
                set(gca,'Fontsize',ft);
                set(gca,'Yscale','log');
                set(gca,'XGrid','on');
                set(gca,'YMinorTick','off');
                set(gca,'linewidth',1);
                plot(mt_tpm,me_tpm-min(Errs(:)),'r-o','linewidth',1.5)
                plot(mt_anls,me_anls-min(Errs(:)),'b->','linewidth',1.5)
                plot(mt_bcd,me_bcd-min(Errs(:)),'c-*','linewidth',1.5)
                plot(mt_ibcd,me_ibcd-min(Errs(:)),'k-p','linewidth',1.5)
                plot(mt_hals,me_hals-min(Errs(:)),'m-d','linewidth',1.5)
                plot(mt_amu,me_amu-min(Errs(:)),'g-<','linewidth',1.5)
                axis([0,timelimit,1.e-4,0.4]);
                box on;
                legend('TPM','BCD','IBCD','ANLS','HALS','AMU','NumColumns',2,'fontsize',12);
                xlabel('CPU Time(s)','fontsize',ft);
                ylabel('$\hat{\rm E}(t)$','Interpreter','latex','fontsize',15);
                set(gca,'linewidth',1);
                title(data(2:end));
            end
        end
end
