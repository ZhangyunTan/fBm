function imgfbf = SimGFBF(Hursts, Poles, ModType,ConvType,rho)
%%        INPUTS
% Hursts: a vector of Hurst parameters, erreur si dim(Hursts>1)
% Poles: a matrix with size (length(Hurst), 2)
% ModType => Modulation Type : '1' for relative frenquencies, '2' for absolute frequencies
% ConvType => Convolution Type : 'full' or 'same', see matlab conv2 function
% rho => dimension of the square input sample matrix = rho*rho 
%%        OUTPUTS
%%
nh = length(Hursts); 
[np1, np2] = size(Poles);

if np2 ~= 2
    error('Sequence -Poles- must ba a matrix with size Np * 2');
end

if nh ~= np1
    error('Sequences Hursts and 2D Poles must have the same size');
end

lxy=rho; 

xfreqs = 0:lxy-1;
yfreqs = 0:lxy-1;
switch ModType;
    case 1 % frequences normalisées
        xfreqs = xfreqs/lxy; 
        yfreqs = yfreqs/lxy;
    case 2 % non-normaliées
       % rien
    case 3 % autre normalisation
        xfreqs = 2*xfreqs/lxy; 
        yfreqs = 2*yfreqs/lxy;
    otherwise
        error('ModulationType must be: 1 or 2');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% WAITBAR
screensize = get(0,'ScreenSize');
PosVect = [1 screensize(4)/2 screensize(3)/6 screensize(4)/12];
h = waitbar(0,'Please wait, computing ...', 'Position',PosVect);
%%%%%%%%%%%%  
% input fbf with pole at 0
fbf = Brownian_field(Hursts(1),lxy);
% samples of interest
imgfbf=fbf(1:lxy,1:lxy);
% modulation, fbf with pole at (Poles(1,1), Poles(1,2))
imgfbf = imgfbf.*( exp(1i*xfreqs'*Poles(1,1))*exp(1i*yfreqs*Poles(1,2)) );
figure
imagesc(abs(imgfbf))
% colormap(gray)
% SaveFigure2Image('imGray1');
% imwrite(abs(imgfbf),'im1.png','png')
% Save2ColorImage(abs(imgfbf),'imRBG1');
waitbar(1/nh,h)
for m=2:nh
    fbf = Brownian_field(Hursts(m),lxy);
    imtemp=fbf(1:lxy,1:lxy);
    imtemp = imtemp.*( exp(1i*xfreqs'*Poles(m,1))*exp(1i*yfreqs*Poles(m,2)) );
    imgfbf = conv2(imgfbf,imtemp, char(ConvType));
    figure
    imagesc(abs(imgfbf))
%     %colormap(gray)
%     SaveFigure2Image(['imGray' num2str(m)]);
%     %imwrite(uint16(abs(imgfbf)),['im' num2str(m) '.png'],'png')
%     Save2ColorImage(abs(imgfbf),['imRBG' num2str(m)]);
waitbar(m/nh,h)
end
close(h)
%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  For comments, send an email to A. Atto: abatt@univ-savoie.fr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%