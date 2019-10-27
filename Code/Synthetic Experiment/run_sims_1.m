clear all
close all
clc

Svec = [4,5,6 7 10 20 50 60 70 80 90];%[4,5,6,10,15,25,40,50,75,100,150];%[150 200 250 500 600 700];%[5,6,10,15,25,40,50];
m = 30;

for S=Svec
    f_Simulation_IRL_1(S,m)
end