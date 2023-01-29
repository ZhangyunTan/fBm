close all;  clear all;  clc; 

%%
N_FBF=2048; 
Hurst = 0.8;  
alpha_attendu=2*Hurst+2;
Wlevel=7; 
Wreg = 'db7'; 

FBF = Brownian_field(Hurst,N_FBF);  % Générer un fBm par méthode Brownian_field.

figure
imagesc(abs(FBF))
colormap(gray)