clear all
close all
clc

% load('data_Ps_5_5_30_samples.mat')
% load('data_Ps_10_10_20_1000_samples.mat')
% load('data_Ps_7_7_30_1000_samples.mat')
load('data_PPP_5_5_30_1000_samples.mat')
Porig = Pcell;
S = 5;
A = 5;
m = 30;
nsvec = [0.001  0.01 0.1 1 10 100]*1e6;
% nsvec = [0.001  ]*1e6;
for nsamp = nsvec
    %     f_Simulation_IRL_L1test_makedata(S,A,m,nsamp);
%     f_Simulation_IRL_L1makeorigdata(S,A,m,nsamp,0.0032);
%     f_Simulation_IRL_L1makeorigdata_pp(S,A,m,nsamp,0.0032);
    tic;f_Simulation_IRL_L1makesampledata(Porig,S,A,m,nsamp);toc
end