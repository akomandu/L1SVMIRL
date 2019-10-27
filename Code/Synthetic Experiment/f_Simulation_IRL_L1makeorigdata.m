function [] = f_Simulation_IRL_L1makeorigdata(S,A,m,nsamp,beta)
% A = 4;
% S = 4;
% nsamp = 500000;
n1 = 50;
n2 =ceil(nsamp/n1);
tol1 = 1e-7;
gamma = 0.1;
% m = 100;
CC = sqrt(S);%sqrt(S);
% RNg = {};
% RhatNg = {};
RMe = {};
RhatMe = {};

lambda = 0.1;

zer = zeros(1,S);
Pcell = {};
PHatcell = {};
parfor kk = 1:m
    Pcol = [];
    [P,Phat] = f_make_data_sp(A,S,n1,n2);
    for a=2:A; Pcol=[Pcol;(P{a}-P{1})]; end
    while(inchull(zer.',Pcol.')||~check_beta(P,A,S,gamma,beta))
        Pcol = [];
        [P,Phat] = f_make_data_sp(A,S,n1,n2);
        for a=2:A; Pcol=[Pcol;(P{a}-P{1})]; end
    end
    Pcell{kk} = P;
    PHatcell{kk} = Phat;
   
end
% for kk = 1:m
%     RdiffNg(kk) = max(abs(RNg{kk} - RhatNg{kk}));
%     RdiffMe(kk) = max(abs(RMe{kk} - RhatMe{kk}));
% end
save(['Data_Files_L1/data_Ps_',num2str(S),'_',num2str(A),'_',num2str(m),'_',num2str(nsamp),'_samples.mat'],'S','A','m','nsamp','Pcell','PHatcell')
end
