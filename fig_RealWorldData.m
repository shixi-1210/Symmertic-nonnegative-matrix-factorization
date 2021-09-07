% Because the experiment time is too long, so we save the result and draw
% the figure to show the performance.
rng(2021);
for fakeloop = 1
    for i =1:4
        if i==1
            load('result_COIL20_edit.mat');
            i_data = 1;
            timelimit= 20;
            data = 'COIL-20';
        elseif i==2
            load('result_ORL_edit.mat');
            i_data =2;
            data = 'ORL';
        elseif i ==3
            i_data = 3;
            load('result_PIE_edit.mat');
            data = 'PIE';
            
        else
            load('result_TDT2_edit.mat');
            i_data = 4;
            data = 'TDT2';
        end
        ft = 13;
        subplot(2,2,i_data); hold on;
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
        title(data);
    end
end