clear all; 

%%
Case = 4;

N_sample = 256; % input sample size rho*rho
N_conv = 9; % number of fbf convolution mixture
eps = 0.3; % small constant

switch Case
    
    case 1
        Hursts = eps + (1-eps)*rand(N_conv,1); % hurst parameters in [eps, 1-eps]\in]0,1[
        Poles = 0.5*pi.*rand(N_conv,2);
        ModulationType = 1;
        ConvolutionType = 'full';
    case 2
        Hursts =0.05 * ones(9,1);
        Poles = pi/8*[1 2 ; 2 1 ; 1 1 ; 2 2 ; 2 2 ; 2 1 ; 1 2 ; 0 0 ; 0 0] ;
        ModulationType = 2;
        ConvolutionType = 'same';
    case 3
        Hursts =[    0.8570    0.9122    0.4172    0.5577    0.8948    0.4697    0.5278    0.4825    0.9381];
        Poles = 0.5*pi.*rand(N_conv,2);
        ModulationType = 2;
        ConvolutionType = 'same';
    case 4
        Hursts =[0.9 0.85];
        Poles = [0 0 ; 0.9817 0.1963] ;
        ModulationType = 2;
        ConvolutionType = 'same';
    otherwise
        error('Set Case = 1 // or Case 2');
end

im = SimGFBF(Hursts, Poles, ModulationType,ConvolutionType,N_sample);

figure
imagesc(abs(im))
colormap(gray)
