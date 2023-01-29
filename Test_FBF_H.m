clear all;  close all; 
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Information: taille d'image d'origine, le parametre Hurst reel, Niveau
% decomposition...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_im = 512; 
Wlevel=7; 
Wreg = 'db45';
Hurst=0.2;
alpha_attendu=2*Hurst+2; 

NbITER=10; 

alpha_diag_8_mc=zeros(1,NbITER); 
alpha_diag_12_mc=zeros(1,NbITER); 
alpha_diag_16_mc=zeros(1,NbITER); 
alpha_diag_20_mc=zeros(1,NbITER); 
alpha_diag_24_mc=zeros(1,NbITER); 
alpha_diag_28_mc=zeros(1,NbITER); 
alpha_diag_32_mc=zeros(1,NbITER); 
alpha_diag_64_mc=zeros(1,NbITER); 

alpha_pola_8_mc=zeros(1,NbITER);
alpha_pola_12_mc=zeros(1,NbITER);
alpha_pola_16_mc=zeros(1,NbITER);
alpha_pola_20_mc=zeros(1,NbITER);
alpha_pola_24_mc=zeros(1,NbITER);
alpha_pola_28_mc=zeros(1,NbITER);
alpha_pola_32_mc=zeros(1,NbITER);
alpha_pola_64_mc=zeros(1,NbITER);

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

[imWPpsd_m,f1,f2]=WPspectrum2D(imag(im_m),Wlevel,Wreg);

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Etape 2: Localisation le pic: estimation les parametres Uq, Vq
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[I,J]=find(imWPpsd_m(:,:)==max(imWPpsd_m(:))); 
Uq_est=(J-1)*pi/2^Wlevel;
Vq_est=(I-1)*pi/2^Wlevel;

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Etape 3: Demodulation: Ramener le pic au point d'originale 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

im_D = im_m.*( exp(-1i*xfreqs'*Uq_est)*exp(-1i*yfreqs*Vq_est) );

%% Spectrum display d�modulation

[imWPpsd_d,f1,f2]=WPspectrum2D(imag(im_D),Wlevel,Wreg);

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Etape 4: Estimer le Parameter H   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
NbPoints=8;
for k=2:NbPoints+1
    for m=k+1:NbPoints+1
        a_Swpdiag=a_Swpdiag+log(imWPpsd_d(m,m)./imWPpsd_d(k,k))/log( sqrt(f1(m+1,m+1)^2 + f2(m+1,m+1)^2) ./ sqrt(f1(k+1,k+1)^2 + f2(k+1,k+1)^2) );
        ncount=ncount+1;
    end;
end;
ncount;
alpha_diag_8=-a_Swpdiag./ncount;
alpha_diag_8_mc(ITER)=alpha_diag_8;

% Methode Polaire 
ncount=0;
a_Swppol=0;
NbPoints=8;
for k=2:NbPoints+1
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
alpha_pola_8_mc(ITER)=alpha_pola_8;

%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimations sur 12 points 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Methode Diagonale 
ncount=0;
a_Swpdiag=0;
NbPoints=12;
for k=2:NbPoints+1
    for m=k+1:NbPoints+1
        a_Swpdiag=a_Swpdiag+log(imWPpsd_d(m,m)./imWPpsd_d(k,k))/log( sqrt(f1(m+1,m+1)^2 + f2(m+1,m+1)^2) ./ sqrt(f1(k+1,k+1)^2 + f2(k+1,k+1)^2) );
        ncount=ncount+1;
    end;
end;
ncount;
alpha_diag_12=-a_Swpdiag./ncount;
alpha_diag_12_mc(ITER)=alpha_diag_12;

% Methode Polaire 
ncount=0;
a_Swppol=0;
NbPoints=12;
for k=2:NbPoints+1
    for m=k+1:NbPoints+1
      if ~isnan(WPpsdRadial(k))
        if ~isnan(WPpsdRadial(m))
            a_Swppol=a_Swppol+log(WPpsdRadialMean(m)./WPpsdRadialMean(k))/log(m./k);
          ncount=ncount+1;
        end
      end
    end;
end; 
alpha_pola_12=-a_Swppol./ncount;
alpha_pola_12_mc(ITER)=alpha_pola_12;

%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimations sur 16 points 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Methode Diagonale 
ncount=0;
a_Swpdiag=0;
NbPoints=16;
for k=2:NbPoints+1
    for m=k+1:NbPoints+1
        a_Swpdiag=a_Swpdiag+log(imWPpsd_d(m,m)./imWPpsd_d(k,k))/log( sqrt(f1(m+1,m+1)^2 + f2(m+1,m+1)^2) ./ sqrt(f1(k+1,k+1)^2 + f2(k+1,k+1)^2) );
        ncount=ncount+1;
    end;
end;
ncount;
alpha_diag_16=-a_Swpdiag./ncount;
alpha_diag_16_mc(ITER)=alpha_diag_16;

% Methode Polaire 
ncount=0;
a_Swppol=0;
NbPoints=16;
for k=2:NbPoints+1
    for m=k+1:NbPoints+1
      if ~isnan(WPpsdRadial(k))
        if ~isnan(WPpsdRadial(m))
            a_Swppol=a_Swppol+log(WPpsdRadialMean(m)./WPpsdRadialMean(k))/log(m./k);
          ncount=ncount+1;
        end
      end
    end;
end; 
alpha_pola_16=-a_Swppol./ncount;
alpha_pola_16_mc(ITER)=alpha_pola_16;
%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimations sur 20 points 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Methode Diagonale 
ncount=0;
a_Swpdiag=0;
NbPoints=20;
for k=2:NbPoints+1
    for m=k+1:NbPoints+1
        a_Swpdiag=a_Swpdiag+log(imWPpsd_d(m,m)./imWPpsd_d(k,k))/log( sqrt(f1(m+1,m+1)^2 + f2(m+1,m+1)^2) ./ sqrt(f1(k+1,k+1)^2 + f2(k+1,k+1)^2) );
        ncount=ncount+1;
    end;
end;
ncount;
alpha_diag_20=-a_Swpdiag./ncount;
alpha_diag_20_mc(ITER)=alpha_diag_20;
% Methode Polaire 
ncount=0;
a_Swppol=0;
NbPoints=20;
for k=2:NbPoints+1
    for m=k+1:NbPoints+1
      if ~isnan(WPpsdRadial(k))
        if ~isnan(WPpsdRadial(m))
            a_Swppol=a_Swppol+log(WPpsdRadialMean(m)./WPpsdRadialMean(k))/log(m./k);
          ncount=ncount+1;
        end
      end
    end;
end; 
alpha_pola_20=-a_Swppol./ncount;
alpha_pola_20_mc(ITER)=alpha_pola_20;
%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimations sur 24 points 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Methode Diagonale 
ncount=0;
a_Swpdiag=0;
NbPoints=24;
for k=2:NbPoints+1
    for m=k+1:NbPoints+1
        a_Swpdiag=a_Swpdiag+log(imWPpsd_d(m,m)./imWPpsd_d(k,k))/log( sqrt(f1(m+1,m+1)^2 + f2(m+1,m+1)^2) ./ sqrt(f1(k+1,k+1)^2 + f2(k+1,k+1)^2) );
        ncount=ncount+1;
    end;
end;
ncount;
alpha_diag_24=-a_Swpdiag./ncount;
alpha_diag_24_mc(ITER)=alpha_diag_24;
% Methode Polaire 
ncount=0;
a_Swppol=0;
NbPoints=24;
for k=2:NbPoints+1
    for m=k+1:NbPoints+1
      if ~isnan(WPpsdRadial(k))
        if ~isnan(WPpsdRadial(m))
            a_Swppol=a_Swppol+log(WPpsdRadialMean(m)./WPpsdRadialMean(k))/log(m./k);
          ncount=ncount+1;
        end
      end
    end;
end; 
alpha_pola_24=-a_Swppol./ncount;
alpha_pola_24_mc(ITER)=alpha_pola_24;
%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimations sur 28 points 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Methode Diagonale 
ncount=0;
a_Swpdiag=0;
NbPoints=28;
for k=2:NbPoints+1
    for m=k+1:NbPoints+1
        a_Swpdiag=a_Swpdiag+log(imWPpsd_d(m,m)./imWPpsd_d(k,k))/log( sqrt(f1(m+1,m+1)^2 + f2(m+1,m+1)^2) ./ sqrt(f1(k+1,k+1)^2 + f2(k+1,k+1)^2) );
        ncount=ncount+1;
    end;
end;
ncount;
alpha_diag_28=-a_Swpdiag./ncount;
alpha_diag_28_mc(ITER)=alpha_diag_28;
% Methode Polaire 
ncount=0;
a_Swppol=0;
NbPoints=28;
for k=2:NbPoints+1
    for m=k+1:NbPoints+1
      if ~isnan(WPpsdRadial(k))
        if ~isnan(WPpsdRadial(m))
            a_Swppol=a_Swppol+log(WPpsdRadialMean(m)./WPpsdRadialMean(k))/log(m./k);
          ncount=ncount+1;
        end
      end
    end;
end; 
alpha_pola_28=-a_Swppol./ncount;
alpha_pola_28_mc(ITER)=alpha_pola_28;
%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimations sur 32 points 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Methode Diagonale 
ncount=0;
a_Swpdiag=0;
NbPoints=32;
for k=2:NbPoints+1
    for m=k+1:NbPoints+1
        a_Swpdiag=a_Swpdiag+log(imWPpsd_d(m,m)./imWPpsd_d(k,k))/log( sqrt(f1(m+1,m+1)^2 + f2(m+1,m+1)^2) ./ sqrt(f1(k+1,k+1)^2 + f2(k+1,k+1)^2) );
        ncount=ncount+1;
    end;
end;
ncount;
alpha_diag_32=-a_Swpdiag./ncount;
alpha_diag_32_mc(ITER)=alpha_diag_32;
% Methode Polaire 
ncount=0;
a_Swppol=0;
NbPoints=32;
for k=2:NbPoints+1
    for m=k+1:NbPoints+1
      if ~isnan(WPpsdRadial(k))
        if ~isnan(WPpsdRadial(m))
            a_Swppol=a_Swppol+log(WPpsdRadialMean(m)./WPpsdRadialMean(k))/log(m./k);
          ncount=ncount+1;
        end
      end
    end;
end; 
alpha_pola_32=-a_Swppol./ncount;
alpha_pola_32_mc(ITER)=alpha_pola_32;
%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimations sur 64 points 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Methode Diagonale 
ncount=0;
a_Swpdiag=0;
NbPoints=64;
for k=2:NbPoints+1
    for m=k+1:NbPoints+1
        a_Swpdiag=a_Swpdiag+log(imWPpsd_d(m,m)./imWPpsd_d(k,k))/log( sqrt(f1(m+1,m+1)^2 + f2(m+1,m+1)^2) ./ sqrt(f1(k+1,k+1)^2 + f2(k+1,k+1)^2) );
        ncount=ncount+1;
    end;
end;
ncount;
alpha_diag_64=-a_Swpdiag./ncount;
alpha_diag_64_mc(ITER)=alpha_diag_64;
% Methode Polaire 
ncount=0;
a_Swppol=0;
NbPoints=64;
for k=2:NbPoints+1
    for m=k+1:NbPoints+1
      if ~isnan(WPpsdRadial(k))
        if ~isnan(WPpsdRadial(m))
            a_Swppol=a_Swppol+log(WPpsdRadialMean(m)./WPpsdRadialMean(k))/log(m./k);
          ncount=ncount+1;
        end
      end
    end;
end; 
alpha_pola_64=-a_Swppol./ncount;
alpha_pola_64_mc(ITER)=alpha_pola_64;
end;
%%
alpha_diag_8_mc_mean=sum(alpha_diag_8_mc)/NbITER;
alpha_diag_12_mc_mean=sum(alpha_diag_12_mc)/NbITER;
alpha_diag_16_mc_mean=sum(alpha_diag_16_mc)/NbITER;
alpha_diag_20_mc_mean=sum(alpha_diag_20_mc)/NbITER;
alpha_diag_24_mc_mean=sum(alpha_diag_24_mc)/NbITER;
alpha_diag_28_mc_mean=sum(alpha_diag_28_mc)/NbITER;
alpha_diag_32_mc_mean=sum(alpha_diag_32_mc)/NbITER;
alpha_diag_64_mc_mean=sum(alpha_diag_64_mc)/NbITER;

alpha_pola_8_mc_mean=sum(alpha_pola_8_mc)/NbITER;
alpha_pola_12_mc_mean=sum(alpha_pola_12_mc)/NbITER;
alpha_pola_16_mc_mean=sum(alpha_pola_16_mc)/NbITER;
alpha_pola_20_mc_mean=sum(alpha_pola_20_mc)/NbITER;
alpha_pola_24_mc_mean=sum(alpha_pola_24_mc)/NbITER;
alpha_pola_28_mc_mean=sum(alpha_pola_28_mc)/NbITER;
alpha_pola_32_mc_mean=sum(alpha_pola_32_mc)/NbITER;
alpha_pola_64_mc_mean=sum(alpha_pola_64_mc)/NbITER;

Biais_diag_8=abs(alpha_diag_8_mc_mean-alpha_attendu); 
Biais_diag_12=abs(alpha_diag_12_mc_mean-alpha_attendu); 
Biais_diag_16=abs(alpha_diag_16_mc_mean-alpha_attendu); 
Biais_diag_20=abs(alpha_diag_20_mc_mean-alpha_attendu); 
Biais_diag_24=abs(alpha_diag_24_mc_mean-alpha_attendu); 
Biais_diag_28=abs(alpha_diag_28_mc_mean-alpha_attendu); 
Biais_diag_32=abs(alpha_diag_32_mc_mean-alpha_attendu); 
Biais_diag_64=abs(alpha_diag_64_mc_mean-alpha_attendu); 

Biais_pola_8=abs(alpha_pola_8_mc_mean-alpha_attendu); 
Biais_pola_12=abs(alpha_pola_12_mc_mean-alpha_attendu); 
Biais_pola_16=abs(alpha_pola_16_mc_mean-alpha_attendu); 
Biais_pola_20=abs(alpha_pola_20_mc_mean-alpha_attendu); 
Biais_pola_24=abs(alpha_pola_24_mc_mean-alpha_attendu); 
Biais_pola_28=abs(alpha_pola_28_mc_mean-alpha_attendu); 
Biais_pola_32=abs(alpha_pola_32_mc_mean-alpha_attendu); 
Biais_pola_64=abs(alpha_pola_64_mc_mean-alpha_attendu); 

Variance_diag_8=sum(alpha_diag_8_mc*alpha_diag_8_mc')/NbITER-alpha_diag_8_mc_mean*alpha_diag_8_mc_mean;
Variance_diag_12=sum(alpha_diag_12_mc*alpha_diag_12_mc')/NbITER-alpha_diag_12_mc_mean*alpha_diag_12_mc_mean;
Variance_diag_16=sum(alpha_diag_16_mc*alpha_diag_16_mc')/NbITER-alpha_diag_16_mc_mean*alpha_diag_16_mc_mean;
Variance_diag_20=sum(alpha_diag_20_mc*alpha_diag_20_mc')/NbITER-alpha_diag_20_mc_mean*alpha_diag_20_mc_mean;
Variance_diag_24=sum(alpha_diag_24_mc*alpha_diag_24_mc')/NbITER-alpha_diag_24_mc_mean*alpha_diag_24_mc_mean;
Variance_diag_28=sum(alpha_diag_28_mc*alpha_diag_28_mc')/NbITER-alpha_diag_28_mc_mean*alpha_diag_28_mc_mean;
Variance_diag_32=sum(alpha_diag_32_mc*alpha_diag_32_mc')/NbITER-alpha_diag_32_mc_mean*alpha_diag_32_mc_mean;
Variance_diag_64=sum(alpha_diag_64_mc*alpha_diag_64_mc')/NbITER-alpha_diag_64_mc_mean*alpha_diag_64_mc_mean;

Variance_pola_8=sum(alpha_pola_8_mc*alpha_pola_8_mc')/NbITER-alpha_pola_8_mc_mean*alpha_pola_8_mc_mean;
Variance_pola_12=sum(alpha_pola_12_mc*alpha_pola_12_mc')/NbITER-alpha_pola_12_mc_mean*alpha_pola_12_mc_mean;
Variance_pola_16=sum(alpha_pola_16_mc*alpha_pola_16_mc')/NbITER-alpha_pola_16_mc_mean*alpha_pola_16_mc_mean;
Variance_pola_20=sum(alpha_pola_20_mc*alpha_pola_20_mc')/NbITER-alpha_pola_20_mc_mean*alpha_pola_20_mc_mean;
Variance_pola_24=sum(alpha_pola_24_mc*alpha_pola_24_mc')/NbITER-alpha_pola_24_mc_mean*alpha_pola_24_mc_mean;
Variance_pola_28=sum(alpha_pola_28_mc*alpha_pola_28_mc')/NbITER-alpha_pola_28_mc_mean*alpha_pola_28_mc_mean;
Variance_pola_32=sum(alpha_pola_32_mc*alpha_pola_32_mc')/NbITER-alpha_pola_32_mc_mean*alpha_pola_32_mc_mean;
Variance_pola_64=sum(alpha_pola_64_mc*alpha_pola_64_mc')/NbITER-alpha_pola_64_mc_mean*alpha_pola_64_mc_mean;

A=[alpha_attendu
alpha_diag_8_mc_mean
alpha_diag_12_mc_mean
alpha_diag_16_mc_mean
alpha_diag_20_mc_mean
alpha_diag_24_mc_mean
alpha_diag_28_mc_mean
alpha_diag_32_mc_mean
alpha_diag_64_mc_mean

alpha_pola_8_mc_mean
alpha_pola_12_mc_mean
alpha_pola_16_mc_mean
alpha_pola_20_mc_mean
alpha_pola_24_mc_mean
alpha_pola_28_mc_mean
alpha_pola_32_mc_mean
alpha_pola_64_mc_mean

Biais_diag_8
Biais_diag_12
Biais_diag_16
Biais_diag_20
Biais_diag_24
Biais_diag_28
Biais_diag_32
Biais_diag_64

Biais_pola_8
Biais_pola_12
Biais_pola_16
Biais_pola_20
Biais_pola_24
Biais_pola_28
Biais_pola_32
Biais_pola_64

Variance_diag_8
Variance_diag_12
Variance_diag_16
Variance_diag_20
Variance_diag_24
Variance_diag_28
Variance_diag_32
Variance_diag_64

Variance_pola_8
Variance_pola_12
Variance_pola_16
Variance_pola_20
Variance_pola_24
Variance_pola_28
Variance_pola_32
Variance_pola_64
]


  

