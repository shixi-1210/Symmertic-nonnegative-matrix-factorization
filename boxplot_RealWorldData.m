% Because the experiment time is too long, so we save the result and draw
% the figure to show the performance.
rng(2021);
for fakeloop = 1
    %绘制箱型图
    for i = 1:2
        figure(i);
        for j =1:2
            k = 2*i+j-2;
            if i == 1 && j==1 %k=1
                name = 'COIL-20';
                load result_COIL20_edit.mat;
            elseif i == 1 && j == 2
                name = 'ORL';
                load('result_ORL_edit.mat');
            elseif i == 2 && j == 1
                name = 'PIE';
                load('result_PIE_edit.mat');
            else
                name = 'TDT2';
                load('result_TDT2_edit.mat');
            end
            EVA = zeros(20,6);
            EVA(:,2) = Clus(k,:,2,1);
            EVA(:,3) = Clus(k,:,3,1);
            EVA(:,4) = Clus(k,:,4,1);
            EVA(:,5) = Clus(k,:,5,1);
            EVA(:,6) = Clus(k,:,6,1);
            EVA(:,1) = Clus(k,:,1,1);
            NAME = {'TPM','ANLS','BCD','IBCD','HALS','AMU'};
            subplot(2,2,2*j-1);
            if 2*j-1 == 1
                subplot('Position',[0.15 0.59 0.25 0.25*1.5]);
            elseif 2*j-1 == 2
                subplot('Position',[0.46 0.59 0.25 0.25*1.5]);
            elseif 2*j-1 == 3
                subplot('Position',[0.15 0.12  0.25 0.25*1.5]);
            elseif 2*j-1 == 4
                subplot('Position',[0.46 0.12  0.25 0.25*1.5]);
            end
            hold on;
            boxplot(EVA,NAME,'color','');
            box on;
            ft = 12;
            set(gca,'Fontsize',ft);
            set(gca,'XGrid','on');
            set(gca,'YMinorTick','off');
            %xlabel('Algorithms','Fontsize',12);
            ylabel('AC','Fontsize',13);
            title(name,'Fontsize',14);
            
            EVA = zeros(20,6);
            EVA(:,1) = Clus(k,:,1,2);
            EVA(:,2) = Clus(k,:,2,2);
            EVA(:,3) = Clus(k,:,3,2);
            EVA(:,4) = Clus(k,:,4,2);
            EVA(:,5) = Clus(k,:,5,2);
            EVA(:,6) = Clus(k,:,6,2);
            subplot(2,2,2*j);
            if 2*j == 1
                subplot('Position',[0.15 0.59 0.25 0.25*1.5]);
            elseif 2*j == 2
                subplot('Position',[0.46 0.59 0.25 0.25*1.5]);
            elseif 2*j == 3
                subplot('Position',[0.15 0.12  0.25 0.25*1.5]);
            elseif 2*j == 4
                subplot('Position',[0.46 0.12  0.25 0.25*1.5]);
            end
            hold on;
            boxplot(EVA,NAME,'color','');
            box on;
            set(gca,'XGrid','on');
            set(gca,'YMinorTick','off');
            set(gca,'Fontsize',ft);
            %xlabel('Algorithms','Fontsize',12);
            ylabel('NMI','Fontsize',13);
            title(name,'Fontsize',14);
            hold on;
        end
    end
end