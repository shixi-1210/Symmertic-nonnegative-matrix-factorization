% This is a figure code of synetic data experiment. Since the synetic world experiment time is too long, we save the result and draw
% the figure to show the performance.
rng(2021);
for fakeloop = 1
    load('Syn_job1');
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