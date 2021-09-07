
clear; close all;
load ('ORL.mat'); 
M = data;
clear data;
truelabel = label;
randn('seed',2018);rand('seed',2018)
X = calcu_similarity_matrix( M );
[m,n] = size(X); r = length(unique(truelabel)); 
U0 = rand(m,r); V0 = U0;
maxiter = 5e2 ; 
lambda = 0.66;  



% --------pivoting ANLS----------
[Up,Vp,ep,tp,diffh,accp] = PivotANLS(X,U0,V0,maxiter,lambda,truelabel);
norm(Up-Vp,'fro')/norm(Vp,'fro')



%---------HALS----------
[Uh,Vh,eh,th,diffh,acch] = HALS_sparse(X,U0,V0,maxiter,lambda,truelabel); %for large sparse matrix
norm(Uh-Vh,'fro')/norm(Vh,'fro')







 % plot 
figure(1)
plot(accp,'r-.','LineWidth',1.5);hold on
plot(acch,'m','LineWidth',1.5);
xlim([0 5e2])
ylim([0.3 1])
% set(gcf,'Position',[100 300 650 300]);
set(gcf,'color','w');
set(gca, 'LineWidth' , 1.5,'FontSize',18);
xlabel('Iteration','FontSize',18);
ylabel('Clustering accuracy','FontSize',18);
legend('SymANLS','SymHALS')


%
figure(2)
plot(tp,accp,'r-.','LineWidth',1.5);hold on
plot(th,acch,'m','LineWidth',1.5);
% set(gcf,'Position',[100 300 650 300]);
set(gcf,'color','w');
set(gca, 'LineWidth' , 1.5,'FontSize',18);
xlim([0 10])
ylim([0.3 1])
xlabel('Time(s)','FontSize',18);
ylabel('Clustering accuracy','FontSize',18);
legend('SymANLS','SymHALS')
