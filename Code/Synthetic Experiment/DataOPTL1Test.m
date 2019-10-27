clear all
close all
clc

% load('data_Ps_5_5_30_samples.mat')
S = 5;
A = 5;
m = 30;
% nsvec = [0.001  ]*1e6;
nsvec = [0.001  0.01 0.1 1 10 100]*1e6;
% nsvec = [  0.01 0.1 1 10 100]*1e6;
for nsamp = nsvec
    %     f_Simulation_IRL_L1test_makedata(S,A,m,nsamp);
    load(['Data_Files_L1/data_Ps2_5_5_30_',num2str(nsamp),'_samples.mat'])
    f_Simulation_IRL_L1Phat(Pcell,PHatcell, S,A,m,nsamp)
%     tic;f_Simulation_IRL_L1makesampledata(Porig,S,A,m,nsamp);toc
end