clear all 
close all
clc

n1 = 50;
n2 = 100;
A = 4;
S = 4;
% gamma = 0.1;
statecell = {};

clear P
for a=1:A; P{a}=randn(S); P{a}=exp(P{a}); P{a}=P{a}./sum(P{a},2); end
for m = 1:n2
rint = randi([1 S], n1+1,1);
states = zeros(n1,3);
x = zeros(n1+1,1);
%% Select initial point
x(1)  = rint(1);
for k = 1:50
%     select action
    temp = simulate(dtmc(P{rint(k+1)}), 1,'X0',full(sparse(1,x(k),1,1,S)));
    x(k+1) = temp(2);
end 
states = [x(1:end-1),rint(2:end), x(2:end)];
statecell{m} = states;
end 

%% Estimate P
PhatCell = estmateP(statecell,S,A);


    
