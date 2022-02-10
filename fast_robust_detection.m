%% Fast and robust pattern detection: Application tospherical bead location in holographic microscopy
% Dylan Brault∗, Loïc Denis∗, Sophie Dixneuf†, Thomas Olivier∗, Nicolas Faure‡, Corinne Fournier∗
% ∗ Univ. Lyon, UJM-Saint-Etienne, CNRS, Institut d Optique Graduate School,Laboratoire Hubert Curien UMR 5516, F-42023, Saint-Etienne, France
% † Bioassays, Microsystems & Optical Engineering Unit, BIOASTER, 69007 Lyon, France
% ‡ bioMerieux, Centre Christophe Mérieux, 38024 Grenoble, France

% (loic.denis@univ-st-etienne.fr, dylan.brault@univ-st-etienne.fr)

% This function computes the likelihood ratio between two hypoethesis :
% - H0 : no object
% - H1 : object at (x,y)
% to perform robust detect of the pattern in an image.

% Inputs :
% f : Robust cost function
% data : Image in which the pattern must be detected
% pattern : Pattern to detect
% k_max : Rank of the low rank approximation used in the computation of the
% robust detection map

%%

function [LR]=fast_robust_detection(f,data,pattern,k_max)
% Correlation function
correl = @(a,b) fftshift(ifft2(fft2(a).*conj(fft2(b))));


% Sampling the possible values of f_s
[A,B] = meshgrid(linspace(min(data(:)),max(data(:)),1000),linspace(min(pattern(:)),max(pattern(:)),1000));
A = A';
B = B';

% Singular value decomposition of F
[U,S,V] = svd(f(A)-f(A-B));
    
LR = zeros(size(data));

% Low rank approximation of the cost function
for k=1:k_max
    LR = LR + correl(reshape(interp1(A(:,1),U(:,k)*sqrt(S(k,k)),data(:),'spline'),size(data)),...
        reshape(interp1(B(1,:),V(:,k)*sqrt(S(k,k)),pattern(:),'spline'),size(data)));
end
LR=fftshift(LR);
end
