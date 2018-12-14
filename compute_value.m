%% 
clear all
close all
a0=0.246e-9;
e=1.602e-19;
t0=2.8*e; % 2.8 ev
n2D=1.2e16; % 1.2 e12 cm-2
h=6.62607e-34;
h_bar=h/(2*pi);
EF=sqrt(9*pi/8)*t0*a0*sqrt(n2D)/e % to get in eV

sf_lim=3*t0*pi/EF/e

R_ext=350;
W=150;
R_mean=(R_ext+(R_ext-W))/2
%delta_B=2*h_bar/e/(271*10^(-9))^2
delta_B=2*h_bar/e/(R_mean*10^(-9))^2