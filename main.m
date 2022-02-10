%% Fast and robust pattern detection: Application tospherical bead location in holographic microscopy
% Dylan Brault∗, Loïc Denis∗, Sophie Dixneuf†, Thomas Olivier∗, Nicolas Faure‡, Corinne Fournier∗
% ∗ Univ. Lyon, UJM-Saint-Etienne, CNRS, Institut d Optique Graduate School,Laboratoire Hubert Curien UMR 5516, F-42023, Saint-Etienne, France
% † Bioassays, Microsystems & Optical Engineering Unit, BIOASTER, 69007 Lyon, France
% ‡ bioMerieux, Centre Christophe Mérieux, 38024 Grenoble, France

% (loic.denis@univ-st-etienne.fr, dylan.brault@univ-st-etienne.fr)
%% Loading image and pattern

% Load the image
data=double(imread('Example/data.png'))/255;

% Load the pattern
pattern=double(imread('Example/pattern.png'))/255;

figure;
subplot(1,2,1);imshow(data,[]);colorbar;title('Data')
subplot(1,2,2);imshow(pattern,[]);colorbar;title('Detection pattern')

%% Robust cost function
% Define the robust cost function
s=mad(data(:),1);

f  = @(x) log(x.^2/s^2+1); % Cauchy loss function
% f = @(x) 1/2*x.^2.*(abs(x)<=s)+s*(abs(x)-1/2*s).*(abs(x)>s); % Huber


k_max=5; % Number of correlation in the approximation

robust_LR=fast_robust_detection(f,data,pattern,k_max);

%% Least-square cost function

correl = @(a,b) ifft2(fft2(a).*conj(fft2(b))); % Correlation function

leastsquare_LR=correl(data,pattern);
figure;

subplot(1,3,1);imshow(fftshift(leastsquare_LR),[]);colormap('jet');colorbar;title('Correlation')
subplot(1,3,2);imshow(fftshift(robust_LR),[]);colormap('jet');colorbar;title('Robust correlation')
subplot(1,3,3);imshow(imread('Example/sample_phase.png'),[]);colorbar;title('Phase of the sample')