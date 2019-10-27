function [] = f_Simulation_IRL_SA(S,A,m)
% A = 4;
% S = 4;
n1 = 50;
n2 =100;
gamma = 0.1;
% m = 100;
CC = sqrt(S);%sqrt(S);
% RNg = {};
% RhatNg = {};
RMe = {};
RhatMe = {};
C = 5;
lambda = 0.1;
cvx_str = make_cvx_str(A);
for kk = 1:m
    [P,Phat] = f_make_data_sp(A,S,n1,n2);

%% My method 
I1 = pinv(eye(S)-gamma*P{1});
clear F
for a=2:A; F{a}=(P{1}-P{a})*I1; end

eval(cvx_str)
RMe{kk} = R;
F2{kk} = P{1};
% hat problem

I1 = pinv(eye(S)-gamma*Phat{1});
clear F
for a=2:A; F{a}=(Phat{1}-Phat{a})*I1; end

eval(cvx_str)
RhatMe{kk} = R;
% save(['Data_Results/data_',num2str(S),'_',num2str(m),'.mat'],'S','m','RMe','RhatMe')
end


% for kk = 1:m
%     RdiffNg(kk) = max(abs(RNg{kk} - RhatNg{kk}));
%     RdiffMe(kk) = max(abs(RMe{kk} - RhatMe{kk}));
% end
save(['Data_Results_SpA/data_',num2str(S),'_',num2str(A),'_',num2str(m),'.mat'],'S','A','m','RMe','RhatMe','F2')

end
