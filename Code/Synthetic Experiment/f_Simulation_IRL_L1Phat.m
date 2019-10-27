function [] = f_Simulation_IRL_L1Phat(Pcell,PHatcell, S,A,m,nsamp)
% A = 4;
% S = 4;
n1 = 50;
n2 =10000;
tol1 = 1e-7;
gamma = 0.1;
% m = 100;
CC = sqrt(S);%sqrt(S);
% RNg = {};
% RhatNg = {};
RMe = {};
RhatMe = {};
C = 5;
lambda = 0.1;
% cvx_str = make_cvx_str_test(A);
zer = zeros(1,S);
for kk = 1:m
    %     Pcol = [];
    %     [P,Phat] = f_make_data_sp(A,S,n1,n2);
    %     for a=2:A; Pcol=[Pcol;P{a}]; end
    %     while(inchull(zer.',Pcol.'))
    %         Pcol = [];
    %         [P,Phat] = f_make_data_sp(A,S,n1,n2);
    %         for a=2:A; Pcol=[Pcol;P{a}]; end
    %     end
    
    %% My method
    P = Pcell{kk};
    I1 = pinv(eye(S)-gamma*P{1});
    
    clear F
    Fcol = [];
    for a=2:A; F{a}=(P{1}-P{a})*I1; end
    for a=2:A; Fcol=[Fcol;F{a}]; end
    %     eval(cvx_str)
    %     RMe{kk} = R;
    %     F2{kk} = P{1};
    FcolOrig{kk} = Fcol;
    
    Phat = PHatcell{kk};
    % hat problem
    I1 = pinv(eye(S)-gamma*Phat{1});
    clear F
    Fcol = [];
    for a=2:A; F{a}=(Phat{1}-Phat{a})*I1; end
    for a=2:A; Fcol=[Fcol;F{a}]; end
    %     eval(cvx_str)
    %     RhatMe{kk} = R;
    FcolHat{kk} = Fcol;
    % save(['Data_Results/data_',num2str(S),'_',num2str(m),'.mat'],'S','m','RMe','RhatMe')
    
    %%
    
    %     I1 = inv(eye(S)-gamma*P{1});
    %     clear F
    %     for a=2:A; F{a}=(P{1}-P{a})*I1; end
    %
    %     cvx_begin
    %     variables R(S) rho
    %     % minimize( 1/2*R.'*R + C * sum((-F{2}-F{3}-F{4}-F{5})*R) )
    %     minimize( lambda*norm(R,1) + sum(max(max(-F{2}*R,-F{3}*R),-F{4}*R)) )
    %     % minimize( 1/2*R.'*R -rho)
    %     subject to
    %     F{2}*R >= zeros(S,1)
    %     F{3}*R >= zeros(S,1)
    %     F{4}*R >= zeros(S,1)
    %
    %     R>= zeros(S,1)-1
    %     R<= zeros(S,1)+1
    %     cvx_end
    %     RNg{kk} = R;
    
    
    I1 = inv(eye(S)-gamma*Phat{1});
    clear F
    for a=2:A; F{a}=(Phat{1}-Phat{a})*I1; end
    
    cvx_begin
    variables R(S) rho
    % minimize( 1/2*R.'*R + C * sum((-F{2}-F{3}-F{4}-F{5})*R) )
    minimize( lambda*norm(R,1) + sum(max(max(max(-F{2}*R,-F{3}*R),-F{4}*R),-F{5}*R)) )
    % minimize( 1/2*R.'*R -rho)
    subject to
    F{2}*R >= zeros(S,1)
    F{3}*R >= zeros(S,1)
    F{4}*R >= zeros(S,1)
    F{5}*R >= zeros(S,1)
    
    R>= zeros(S,1)-1
    R<= zeros(S,1)+1
    cvx_end
    RhatNg{kk} = R;
    
    cvx_begin
    variables R(S)
    % minimize( 1/2*R.'*R + C * sum((-F{2}-F{3}-F{4}-F{5})*R) )
    minimize( norm(R,1))
    % minimize( 1/2*R.'*R -rho)
    subject to
    F{2}*R >= ones(S,1)
    F{3}*R >= ones(S,1)
    F{4}*R >= ones(S,1)
    F{5}*R >= ones(S,1)
    
    
    cvx_end
    RhatSVM{kk} = [R];
end
% for kk = 1:m
%     RdiffNg(kk) = max(abs(RNg{kk} - RhatNg{kk}));
%     RdiffMe(kk) = max(abs(RMe{kk} - RhatMe{kk}));
% end
save(['Data_Files_L1/results_L1_',num2str(S),'_',num2str(A),'_',num2str(m),num2str(nsamp),'.mat'],'S','A','m','RhatSVM','RhatNg','FcolOrig','FcolHat')

end
