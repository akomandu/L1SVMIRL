function [] = f_Simulation_IRL_L1makesampledata(Porig,S,A,m,nsamp)
% A = 4;
% S = 4;
% nsamp = 500000;
n1 = 100;
n2 =ceil(nsamp/n1);
tol1 = 1e-7;

% RNg = {};
% RhatNg = {};





Pcell = {};
PHatcell = {};
parfor kk = 1:m
    [Phat] = f_make_sample_p(Porig{kk},A,S,n1,n2);
    Pcell{kk} = Porig{kk};
    PHatcell{kk} = Phat;
   
end
% for kk = 1:m
%     RdiffNg(kk) = max(abs(RNg{kk} - RhatNg{kk}));
%     RdiffMe(kk) = max(abs(RMe{kk} - RhatMe{kk}));
% end
save(['Data_Files_L1/data_Ps2_',num2str(S),'_',num2str(A),'_',num2str(m),'_',num2str(nsamp),'_samples.mat'],'S','A','m','nsamp','Pcell','PHatcell')
end
