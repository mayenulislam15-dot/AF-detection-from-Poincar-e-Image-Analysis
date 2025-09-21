%% 2D DWT
originalImage = resized_image;
wavelet = 'db10';  % Daubechies wavelet with 1 vanishing moment (change as needed)
levels = 2;        % Number of decomposition levels
[C, S] = wavedec2(originalImage, levels,wavelet);

%level1 decomposition
[H1,V1,D1] = detcoef2('all',C,S,1);
A1 = appcoef2(C,S,'db10',1);
% %level 2 decomposition
[H2,V2,D2] = detcoef2('all',C,S,2);
A2 = appcoef2(C,S,'db10',2);
%Spectral Energy features
EA1 = sum(A1(:).^2);
EA2 = sum(A2(:).^2);
EH1 = sum(H1(:).^2);
EH2 = sum(H2(:).^2);
EV1 = sum(V1(:).^2);
EV2 = sum(V2(:).^2);
ED1 = sum(D1(:).^2);
ED2 = sum(D2(:).^2);

%Entropy Features
Entr_A1 = calculateEntropy(A1);
Entr_A2 = calculateEntropy(A2);
Entr_H1 = calculateEntropy(H1);
Entr_H2 = calculateEntropy(H2);
Entr_V1 = calculateEntropy(V1);
Entr_V2 = calculateEntropy(V2);
Entr_D1 = calculateEntropy(D1);
Entr_D2 = calculateEntropy(D2);

%Normalized entropy
Nor_Entr_A1 = calculateEntropy(mat2gray(A1));
Nor_Entr_A2 = calculateEntropy(mat2gray(A2));
Nor_Entr_H1 = calculateEntropy(mat2gray(H1));
Nor_Entr_H2 = calculateEntropy(mat2gray(H2));
Nor_Entr_V1 = calculateEntropy(mat2gray(V1));
Nor_Entr_V2 = calculateEntropy(mat2gray(V2));
Nor_Entr_D1 = calculateEntropy(mat2gray(D1));
Nor_Entr_D2 = calculateEntropy(mat2gray(D2));

%% 2D EWT 
params.globtrend = 'none';
params.degree=5; 
params.reg = 'gaussian';
params.lengthFilter = 10;
params.sigmaFilter = 1.5;
params.detect = 'scalespace';
params.typeDetect='otsu'; 
params.N = 3; 
params.completion=0;
params.log=0;
Bound=1;   
Comp=1;    
Rec=1;     
[ewtc,mfbR,mfbC,BR,BC]=EWT2D_Tensor(rgb2gray(resized_image),params);
ewtc1 = ewtc{1,1};
[md_ewt, nd_ewt] = size(ewtc);
ewtcf = ewtc{max(md_ewt), max(nd_ewt)};
%%2D EWT Features
E_ewt1 = sum(ewtc1(:).^2); % Spectral energy of ewtc1
E_ewtf = sum(ewtcf(:).^2); % Spectral energy of ewtcf
Entr_ewt1 = calculateEntropy(ewtc1);
Entr_ewtf = calculateEntropy(ewtcf);
Entr_ewt1_nor = calculateEntropy(mat2gray(ewtc1));
Entr_ewtf_nor = calculateEntropy(mat2gray(ewtcf));

%% Curvelet transform
c_img= rgb2gray(resized_image);
Curvelet = fdct_wrapping(double(c_img), 1, 2, 3);
cA = Curvelet{1,1}{1,1}; 
cAd = 255.*(cA./max(cA)); 
%%Curvelet features
Energy_Appcoeff = sum(cA(:).^2); 
Entropy_Appcoeff = calculateEntropy(cA);
Entropy_nor_Appcoeff = calculateEntropy(mat2gray(cA));
Energy_original = sum(Curvelet{1,3}{1,1}(:).^2); 
Entropy_original = calculateEntropy(Curvelet{1,3}{1,1});
Entropy_nor_original = calculateEntropy(mat2gray(Curvelet{1,3}{1,1}));
for k_curvelet = 1:16
    detail_coeff=Curvelet{1,2}{1,k_curvelet};
    Energy_detailcoeff(k_curvelet) = sum(Curvelet{1,2}{1,k_curvelet}(:).^2);
    Entropy_detailcoeff(k_curvelet) =  calculateEntropy(Curvelet{1,2}{1,k_curvelet});
    Entropy_nor_detailcoeff(k_curvelet) = calculateEntropy(mat2gray(Curvelet{1,2}{1,k_curvelet}));
end


