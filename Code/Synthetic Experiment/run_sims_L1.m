clear all
close all
clc

Avec = [4];%[4,5,6,10,15,25,40,50,75,100,150];%[150 200 250 500 600 700];%[5,6,10,15,25,40,50];
Svec = [5];%,7,10,12,15];
m = 30;
for S = Svec
for A=Avec
    f_Simulation_IRL_L1test(S,A,m)
end
end