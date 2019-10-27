function [P,Phat] = f_make_data_sp(A,S,n1,n2)
K = 1;
pr1 = 0.2;
% n1 = 50;
% n2 = 100;
% A = 4;
% S = 4;
% gamma = 0.1;
statecell = {};

clear P
for a=1:A; Pt{a}=randn(S); Pt{a}=exp(Pt{a}); Pt{a}=Pt{a}./sum(Pt{a},2);
    k1 = [];
    k2 = [];
    for b =  1:S
        k1 = [k1;b*ones(K,1)];
        k2 = [k2;randperm(S,K).'];
    end
    Q{a}=randn(S,K); Q{a}=exp(Q{a}); Q{a}=[Q{a}./sum(Q{a},2)].';
    QQ{a} = sparse(k1,k2,Q{a}(:),S,S);
    
    P{a} = pr1*Pt{a} + (1-pr1)*QQ{a};
    
end

for m = 1:n2
rint = [randi([1 S]);randi([1 A], n1,1)];
states = zeros(n1,3);
x = zeros(n1+1,1);
%% Select initial point
x(1)  = rint(1);
for k = 1:n1
%     select action
    temp = simulate(dtmc(P{rint(k+1)}), 1,'X0',full(sparse(1,x(k),1,1,S)));
    x(k+1) = temp(2);
end 
states = [x(1:end-1),rint(2:end), x(2:end)];
statecell{m} = states;
end 

%% Estimate P
Phat = estmateP(statecell,S,A);


    
