clear all;  close all;clc; 
%%

%   im_o = imread('D:\Documents\MET\MET_Images_11_09_2014\0907241-001-014C400-140.bmp');
  im_o = imread('0907241-001-03CROP2.bmp');
%   im_o = imread('D:\Documents\MET\MET_Images_11_09_2014\0907241-001-010C290-360.bmp');
%  im_o = imread('D:\Documents\MET\MET_Images_11_09_2014\0907241-001-013C280-60.bmp');
im_o = double(im_o(:,:,1));
im55 = im_o - mean(mean(im_o));
im=im55;

% Interpolation 
%im = imresize(im55,[128,128]);

%im = imrotate(im55,90);

dim_I = size(im);
dim = dim_I(1);

% Module FFT
Puiss_Mod = 1.5;
ModFFT_im = fftshift(abs(fft2(im)));

figure
subplot(1,2,1)
imagesc(im)
axis image;
colormap(gray)
title('Image HRTEM')
subplot(1,2,2)
imagesc(ModFFT_im.^Puiss_Mod)
axis image;
title('Periodogram')
%%
J0=4;
Wreg='db7';

%%%%%%%%%%%%%%%%%%%%%%%%
% WPpsd_im_o

 [WPpsd_im,w1,w2,alpha_im]=WP2D_Estim_Spectrum_IsotAlpha(im,J0,Wreg, 2);
%  [WPpsd_im,w1,w2,alpha_im]=WP2D_Estim_Spectrum_IsotAlpha_1(im,J0,Wreg);

figureFS = figure('Color',[1 1 1]);
axes3 = axes('Visible','off','Parent',figureFS,'CLim',[0 1]);
surfl(w1,w2,WPpsd_im);
set(gcf,'Color','w')
hold on
stem3(w1,w2,WPpsd_im,'k','MarkerSize',3);
view(axes3,[-225 20]);
title('WPpsd of original image - Full --> Support = [0,\pi]\times [0,\pi]')

if (J0>4)
    Ns1 = round(2^J0/2); % number of x-axis samples : 1:Ns1 ---> [0 pi/2]
    Ns2 = round(2^J0/2); % number of y-axis samples : 1:Ns2 ---> [0 pi/2]
    %%%% ZOOM ROTATE /
    %
    figureZR = figure('Color',[1 1 1]);
    axes2 = axes('Visible','off','Parent',figureZR,'CLim',[0 1]);
    surfl(w1(1:Ns1,1:Ns2),w2(1:Ns1,1:Ns2),WPpsd_im(1:Ns1,1:Ns2),'Parent',axes2);
    grid(axes2,'on');
    hold(axes2,'all');
    %set(gcf,'Color','w')
    view(axes2,[-225 20]);
    title('WP Spectrum démodulation over [0,\pi/2]\times [0,\pi/2]')
    %
    %%%% ZOOM PROFILE /
    figureZP = figure('Color',[1 1 1]);
    axes1 = axes('Visible','off','Parent',figureZP,'CLim',[0 1]);
    surfl(w1(1:Ns1,1:Ns2),w2(1:Ns1,1:Ns2),WPpsd_im(1:Ns1,1:Ns2),'Parent',axes1);
    grid(axes1,'on');
    hold(axes1,'all');
    %set(gcf,'Color','w')
    view(axes1,[-45-180 0]);
    title('WP Spectrum démodulation over [0,\pi/2]\times [0,\pi/2] - Profile view')
end
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Image fBmRemoval

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MyTree=wpdec2(im,J0,Wreg);

PoleProcessing = 2; % CHOIX 1 : BLANCHIMENT, Choix 2, 3, etc: SUPPRESSION /f^alpha
TreeTemp = MyTree;
switch PoleProcessing
    case 1 % Blanchiment
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for k=0:4^J0-1
            imgcoeff=wpcoef(TreeTemp,[J0,k]);
            gamma_jn=std2(imgcoeff);
            imgcoeff=imgcoeff./gamma_jn;
            TreeTemp=write(TreeTemp,'cfs',[J0 k],imgcoeff);
        end
    otherwise % Suppression du pole en 0
        wisot = sqrt(w1.^2 + w2.^2).^(alpha_im/2);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 1) Frequency 2^Wlevel points
        Wlevel = J0;
        p=0:1:2^Wlevel-1;
        n = bitxor(p, fix(p/2));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tab = [1 2; 3 4]; % tab-1 is the quaternary basis for \mu = 0, 1, 2, 3
        for u= 1:Wlevel-1
            tab = [tab tab+4^u];
        end
        for v=1:Wlevel-1
            tab = [tab; tab+2*4^v];
        end
        % tab contains the capital frequency indices ordered with respect to the
        % rectangular grid with size p \times p
        for p1=1:3 % 2^Wlevel
            for p2=1:3 %2^Wlevel
                imgcoeff=wpcoef(TreeTemp,[Wlevel -1+tab(n(p1)+1,n(p2)+1)]);
                
                imgcoeff=imgcoeff.*( wisot(p1,p2) );
                TreeTemp=write(TreeTemp,'cfs',[Wlevel -1+tab(n(p1)+1,n(p2)+1)],imgcoeff);
            end
        end
end
%%%%%%%% RECONSTRUCTION %%%%%%%%%%%%%
I_fBmRemoval=wprec2(TreeTemp);
I_fBmRemoval = double(I_fBmRemoval(:,:,1));
%%
% Module I_fBmRemoval
ModFFT_fBmRemoval = fftshift(abs(fft2(I_fBmRemoval)));

% figure
% subplot(1,2,1)
% imagesc((I_fBmRemoval-min(min(I_fBmRemoval)))/(max(max(I_fBmRemoval))-min(min(I_fBmRemoval)))); 
% axis image; 
% colormap(gray)
% title('HRTEM with fBmRemoval')
% subplot(1,2,2)
% imagesc(ModFFT_fBmRemoval.^Puiss_Mod)
% axis image;
% title('Periodogram')

% WPpsd_fBmRemoval
 [WPpsd_fBmRemoval,w1,w2,alpha_fBmRemoval]=WP2D_Estim_Spectrum_IsotAlpha(I_fBmRemoval,J0,Wreg,10);
%  [WPpsd_fBmRemoval,w1,w2,alpha_fBmRemoval]=WP2D_Estim_Spectrum_IsotAlpha_1(I_fBmRemoval,J0,Wreg);

figureFS = figure('Color',[1 1 1]);
axes3 = axes('Visible','off','Parent',figureFS,'CLim',[0 1]);
surfl(w1,w2,WPpsd_fBmRemoval);
set(gcf,'Color','w')
hold on
stem3(w1,w2,WPpsd_fBmRemoval,'k','MarkerSize',3);
view(axes3,[-225 20]);
title('WPpsd_fBmRemoval - Full --> Support = [0,\pi]\times [0,\pi]')

if (J0>4)
    Ns1 = round(2^J0/2); % number of x-axis samples : 1:Ns1 ---> [0 pi/2]
    Ns2 = round(2^J0/2); % number of y-axis samples : 1:Ns2 ---> [0 pi/2]
    %%%% ZOOM ROTATE /
    %
    figureZR = figure('Color',[1 1 1]);
    axes2 = axes('Visible','off','Parent',figureZR,'CLim',[0 1]);
    surfl(w1(1:Ns1,1:Ns2),w2(1:Ns1,1:Ns2),WPpsd_fBmRemoval(1:Ns1,1:Ns2),'Parent',axes2);
    grid(axes2,'on');
    hold(axes2,'all');
    %set(gcf,'Color','w')
    view(axes2,[-225 20]);
    title('WP Spectrum démodulation over [0,\pi/2]\times [0,\pi/2]')
    %
    %%%% ZOOM PROFILE /
    figureZP = figure('Color',[1 1 1]);
    axes1 = axes('Visible','off','Parent',figureZP,'CLim',[0 1]);
    surfl(w1(1:Ns1,1:Ns2),w2(1:Ns1,1:Ns2),WPpsd_fBmRemoval(1:Ns1,1:Ns2),'Parent',axes1);
    grid(axes1,'on');
    hold(axes1,'all');
    %set(gcf,'Color','w')
    view(axes1,[-45-180 0]);
    title('WP Spectrum démodulation over [0,\pi/2]\times [0,\pi/2] - Profile view')
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Etape 2: Localisation un pic aprtir de l'image I_fBmRemoval
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[I,J]=find(WPpsd_fBmRemoval(:,:)==max(WPpsd_fBmRemoval(:))); 
Uq_est=(J-1)*pi/2^Wlevel
Vq_est=(I-1)*pi/2^Wlevel

% Uq_est=pi/16;
% Vq_est=pi/32;

xfreqs = 0:dim-1;
yfreqs = 0:dim-1;
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Etape 3: Demodulation: Ramener le pic au point d'originale 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
im_D = I_fBmRemoval.*( exp(-1i*xfreqs'*Uq_est)*exp(-1i*yfreqs*Vq_est) );
[imWPpsd_d,f1,f2]=WPspectrum2D(im_D,J0,Wreg);

%% Transformee polaire

[M,N] = size(imWPpsd_d);
imW=zeros(2*M-1,2*N-1);
for j=1:N
  for i=1:M
    imW(i,j)=imWPpsd_d(M-i+1,N-j+1);
  end
end
for j=1:N
  for i=1:M
    imW(i,N+j-1)=imWPpsd_d(M-i+1,j);
  end
end
for j=1:N
  for i=1:M
    imW(M+i-1,N+j-1)=imWPpsd_d(i,j);
  end
end
for j=1:N
  for i=1:M
    imW(M+i-1,j)=imWPpsd_d(i,N-j+1);
  end
end
clear WPpsdRadialMean;

%calcul transformee polaire avec 360 pas angulaires

[WPpsdRadial]=polartrans(imW,N,360);
WPpsdRadialMean = zeros(1,N);
for m=1:N
     WPpsdRadialMean(m) = mean(WPpsdRadial(m,:));
end

%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimations sur 8 points 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Methode Diagonale 
ncount=0;
a_Swpdiag=0;
NbPoints=3;
for k=1:NbPoints+1
    for m=k+1:NbPoints+1
        a_Swpdiag=a_Swpdiag+log(imWPpsd_d(m,m)./imWPpsd_d(k,k))/log( sqrt(f1(m+1,m+1)^2 + f2(m+1,m+1)^2) ./ sqrt(f1(k+1,k+1)^2 + f2(k+1,k+1)^2) );
        ncount=ncount+1;
    end;
end;
ncount;
alpha_diag_2=-a_Swpdiag./ncount;

% Methode Polaire 
ncount=0;
a_Swppol=0;
NbPoints=8;
for k=1:NbPoints+1
    for m=k+1:NbPoints+1
      if ~isnan(WPpsdRadial(k))
        if ~isnan(WPpsdRadial(m))
            a_Swppol=a_Swppol+log(WPpsdRadialMean(m)./WPpsdRadialMean(k))/log(m./k);
          ncount=ncount+1;
        end
      end
    end;
end; 
alpha_pola_8=-a_Swppol./ncount;

A= [
    alpha_im 
    alpha_fBmRemoval 
    alpha_diag_2 
    alpha_pola_8
    ]

% if alpha_im >4
%     alpha_im = 3.99;
% end
% 
% if alpha_fBmRemoval < 0
%     alpha_fBmRemoval = -alpha_fBmRemoval;
% end
%% Synthese image

N_sample = 128; % input sample size rho*rho
N_conv = 9; % number of fbf convolution mixture
eps = 0.3; % small constant

Hursts =[(alpha_im-2)/2 (alpha_fBmRemoval-2)/2];
        Poles = [0 0 ; Uq_est Vq_est] ;
        ModulationType = 2;
        ConvolutionType = 'same';
im_syn = SimGFBF(Hursts, Poles, ModulationType,ConvolutionType,N_sample);
%
figure
imagesc(abs(im_syn))
colormap(gray)

% 
% 
%  [WPpsd_imsyn,w1,w2,alpha_im]=WP2D_Estim_Spectrum_IsotAlpha(im_syn,J0,Wreg, 2);
% %  [WPpsd_im,w1,w2,alpha_im]=WP2D_Estim_Spectrum_IsotAlpha_1(im,J0,Wreg);
% 
% figureFS = figure('Color',[1 1 1]);
% axes3 = axes('Visible','off','Parent',figureFS,'CLim',[0 1]);
% surfl(w1,w2,WPpsd_imsyn);
% set(gcf,'Color','w')
% hold on
% stem3(w1,w2,WPpsd_imsyn,'k','MarkerSize',3);
% view(axes3,[-225 20]);
% title('WPpsd_imsyn- Full --> Support = [0,\pi]\times [0,\pi]')
