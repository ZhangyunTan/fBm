function [WPpsd,w1,w2,alpha]=WP2D_Estim_Spectrum_IsotAlpha(img,Wlevel,Wreg,NbPoints)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computation of the wavelet packet spectrum, see paper
%  ''2D Wavelet Packet Spectrum for Image Analysis'' by A. Atto,
%               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REQUIRES:
% matlab 'wavelet toolbox'
%        INPUTS
% img: input grayscale image
% Wlevel:
%    Wavelet Packets decomposition level, command 'doc wpdec2' for details.
% Wreg:
%    Wavelet name, Example: 'sym8', command 'doc waveinfo' for details.
%        OUTPUTS
% WPpsd_var: wavelet packet spectrum (given over meshgrid w1 and w2).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1) Frequency 2^Wlevel points
p=0:1:2^Wlevel-1;
[w1,w2] = meshgrid(p*pi/2^Wlevel); % frequency grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2)   Compute the capital (quaternary based) frequency indices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reorder frequency points with respect to the Gray code permutation
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Zero-mean textures and random fields are considered (see paper)
img = double(img);
mean2(img);
img = img - mean2(img); %
mean2(img);
WPpsd=zeros(2^Wlevel);
%%
% 3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Wavelet Packet Spectrum at Wlevel with Wreg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
WTree=wpdec2(img,Wlevel,Wreg);
for p1=1:2^Wlevel
    for p2=1:2^Wlevel
        Yjk=wpcoef(WTree,[Wlevel -1+tab(n(p1)+1,n(p2)+1)]);
        WPpsd(p1,p2)=var(Yjk(:));
    end
end

%% Diagonal alpha Estimation for Isotropic 1/f^alpha field
%
%
WP_psd1 = diag(WPpsd);
ncount = 0;
a_Swp1=0;
%NbPoints=4;
for k=1:NbPoints
    for m=k+1:NbPoints
        a_Swp1=a_Swp1+log(WP_psd1(m)./WP_psd1(k))/log((m)/(k));
        ncount=ncount+1 ;
    end;
end;
alpha=-a_Swp1./ncount;
%alpha=sum(a_Swp1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright: Abdourrahmane.Atto@univ-savoie.fr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%