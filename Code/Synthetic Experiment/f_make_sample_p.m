function [Phat] = f_make_sample_p(P,A,S,n1,n2)

% n1 = 50;
% n2 = 100;
% A = 4;
% S = 4;
% gamma = 0.1;
statecell = {};
%% Make large P matrix
Plarge = [];
for ka = 1:A
    Plarge = [Plarge;(repmat(P{ka},1,A))/A];
end

%%

for m = 1:n2
rint = [randi([1 S*A])];
states = zeros(n1,3);
% x = zeros(n1+1,1);
%% Select initial point
x(1)  = rint(1);

%     select action
    x = simulate(dtmc(Plarge), n1,'X0',full(sparse(1,x(1),1,1,S*A)));
    
states = [mod(x(1:end-1),S),ceil(x(1:end-1)/S), mod(x(2:end),S)];
states(states == 0) = S;
statecell{m} = states;
end 

%% Estimate P
Phat = estmateP(statecell,S,A);


    
