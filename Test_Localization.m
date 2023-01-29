clear all; close all;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Information: taille d'image d'origine, le parametre Hurst reel, Niveau
% decomposition...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N_im = 512; 
Wlevel=7; 
Wreg = 'db45';
Hurst = 0.4; 
alpha_attendu=2*Hurst+2; 

NbITER=10; 

Uq_r_mc = zeros(1,NbITER);
Vq_r_mc = zeros(1,NbITER);

Uq_i_mc = zeros(1,NbITER);
Vq_i_mc = zeros(1,NbITER);

for ITER=1:NbITER 
im_o = Brownian_field(Hurst,N_im);
im_o = double(im_o);
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Etape 1: Modulation, FBF avec un pic � zero deplacer au point (Uq,Vq)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Uq = pi/4;
Vq = pi/8;
xfreqs = 0:N_im-1;
yfreqs = 0:N_im-1;    

im_m = im_o.*( exp(1i*xfreqs'*Uq)*exp(1i*yfreqs*Vq) );
%% Spectrum FBF Modulation
[imWPpsd_real,f1,f2]=WPspectrum2D(real(im_m),Wlevel,Wreg);
[imWPpsd_imag,f1,f2]=WPspectrum2D(imag(im_m),Wlevel,Wreg);
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Etape 2: Localisation le pic: estimation les parametres Uq, Vq
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[I_r,J_r]=find(imWPpsd_real(:,:)==max(imWPpsd_real(:))); 
Uq_r=(J_r-1)*pi/2^Wlevel;
Vq_r=(I_r-1)*pi/2^Wlevel;
Uq_r_mc(ITER)=Uq_r;
Vq_r_mc(ITER)=Vq_r;


[I_i,J_i]=find(imWPpsd_imag(:,:)==max(imWPpsd_imag(:))); 
Uq_i=(J_i-1)*pi/2^Wlevel;
Vq_i=(I_i-1)*pi/2^Wlevel;
Uq_i_mc(ITER)=Uq_r;
Vq_i_mc(ITER)=Vq_r;

Uq_r_m=sum(Uq_r_mc)/NbITER;
Vq_r_m=sum(Vq_r_mc)/NbITER;

Uq_i_m=sum(Uq_i_mc)/NbITER;
Vq_i_m=sum(Vq_i_mc)/NbITER;


Uq_r_v=var(Uq_r_mc);
Vq_r_v=var(Vq_r_mc);

Uq_i_v=var(Uq_i_mc);
Vq_i_v=var(Vq_i_mc);
end; 

A=[ 
    Uq Vq
    Uq_r_m Vq_r_m
    Uq_i_m Vq_i_m
    Uq_r_v Vq_r_v
    Uq_i_v Vq_i_v
    ]

