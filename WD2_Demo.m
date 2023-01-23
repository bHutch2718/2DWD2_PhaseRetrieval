%% PHASE RETRIEVAL BY 2DWD2: 2D Wigner Distribution Deconvolution

% This script provides two demonstations of 2D phase retrieval by Wigner
% Distribution Deconvolution.  The algorithm is discussed in the paper:
%
% "Fast and Accurate 2D Phase Retrieval from Incomplete, Windowed Fourier
% Measurements"
% Maya Hamka, Brandon Hutchison, Donald Liveoak
% (link incoming...)
%
% This work was supported in part by the National Science Foundation
% (DMS-2012238), and the Dept. of Mathematics & Statistics, University of
% Michigan - Dearborn.
% Mentors: Profs. Aditya Viswanathan and Yulia Hristova

clear;clc;close;
rng(2718);

%% Query for Type of Signal

fprintf( ' --------------------------------------------------- \n' );
fprintf( ' --- 2DWD2: 2D Wigner Distribution Deconvolution --- \n' );
fprintf( ' --------------------------------------------------- \n\n' );

fprintf( ' This demonstration recovers two complex signal types:\n\n');
fprintf( '     a) Randomly Generated \n');
fprintf( '     b) Deterministic with Test Images \n\n');

prompt = ' Enter signal type: (a or b) ';
type = input(prompt,"s");

% Default to random sig. for bad user input
if (~strcmp(type,'a')) && (~strcmp(type,'b'))
    type = 'a';
end

%% Problem Parameters, Test Signal and Mask Generation

switch type
    case 'a'  % Random Complex Signal
        fprintf( '\n --- Generating measurements for RANDOM SIGNAL --- \n\n');
        % Size of Test Signal and SNR
        d1 = 228;
        d2 = 204;
        SNR = 30;
        
        % Number of Samples
        K1 = 19;
        K2 = 17;
        K1_of_em = 0:d1/K1:d1-1;
        K2_of_em = 0:d2/K2:d2-1;

        % Support of Mask
        delta1 = 10;
        delta2 =9;
        kappa1 = K1 - delta1 + 1;
        kappa2 = K2 - delta2 + 1;
        
        % Generate Signal
        X = randn(d1, d2) + 1i*randn(d1,d2);

    case 'b'  % Deterministic Test Images
        fprintf( '\n --- Generating measurements for TEST IMAGE --- \n\n');
        % Size of Test Signal and SNR
        d1 = 253;
        d2 = 200;
        SNR = 60;
        
        % Number of Samples
        K1 = 11;
        K2 = 5;
        K1_of_em = 0:d1/K1:d1-1;
        K2_of_em = 0:d2/K2:d2-1;

        % Support of Mask
        delta1 = (K1+1)/2;
        delta2 = (K2+1)/2;
        kappa1 = delta1;
        kappa2 = delta2;

        % Set magnitude to image of Fourier
        magImage = imread('image_magFourier.jpg');
        magImage = rgb2gray(magImage);
        magImage = imresize(magImage,[d1 d2]);
        magImage = im2double(magImage);

        % Set phase to image of Rodenburg
        phaseImage = imread('image_phaseRodenburg.png');
        phaseImage = rgb2gray(phaseImage);
        phaseImage = imresize(phaseImage,[d1 d2]);
        phaseImage = im2double(phaseImage);
        
        % Generate Signal
        X = abs(magImage) .* exp(1i.*phaseImage);
end

vectLength = (4*kappa1*kappa2-2*kappa1-2*kappa2+1)*d1*d2;

% Generate Random Mask
M = zeros(d1,d2);
M(1:delta1, 1:delta2) = randn(delta1, delta2) + 1i*randn(delta1,delta2);

% Print Problem Parameters
fprintf( ' Problem size: d1 x d2 = %d x %d \n', d1, d2 );
fprintf( ' Window size: %d x %d \n', delta1, delta2 );
                                 
fprintf( ' No. of Fourier modes: K1 x K2 = %d x %d = %d \n', K1, K2, K1*K2 );
fprintf( ' Total no. of measurements: %d(d1 x d2) = %d \n\n', K1*K2, d1*d2*K1*K2 );

fprintf( ' Signal-to-Noise Ratio: %d dB \n\n', SNR );

%%  Gernerate Measurements and Add Noise

% Measurement Tensor
Y = zeros(K1,K2,d1,d2);

for L1 = 0:d1-1
    for L2 = 0:d2-1
        fullSlice = abs(fft2( X .* circshift(M, [L1,L2]) )).^2;
        Y(:,:,L1+1,L2+1) = fullSlice(K1_of_em +1, K2_of_em +1);
    end
end
clear fullSlice;

% Noise
noiseAmp = norm( Y(:) ) / sqrt( K1*K2*d1*d2*10^(SNR/10) );
N = noiseAmp * randn(K1,K2,d1,d2);
Y = Y + N;
clear N;

%% Mask Pre-Computations

maskNumber = 0;
maskComps = zeros(d1, d2, 4*kappa1*kappa2-2*kappa1-2*kappa2+1);
baseIndex = reshape(1:d1*d2,d1,d2);
sparseIndices = zeros(vectLength,2);

for k1 = 0:kappa1-1
    for k2 = 0:kappa2-1
        maskNumber = maskNumber + 1;
        M_0 = rev2(fft2( M.*circshift( conj(M), -[k1 k2]) ));
        maskComps(:,:, maskNumber) = M_0;
        shiftIndex = circshift(baseIndex, -[k1 k2]);
        sparseIndices((maskNumber-1)*d1*d2+1:maskNumber*d1*d2, :) = [baseIndex(:),shiftIndex(:)]; 
    end
end

for k1 = delta1:K1-1
    for k2 = delta2:K2-1
        maskNumber = maskNumber + 1;
        M_0 = rev2(fft2( M.*circshift( conj(M), [K1-k1 K2-k2]) ));
        maskComps(:,:, maskNumber) = M_0;
        shiftIndex = circshift(baseIndex, [K1-k1 K2-k2]);
        sparseIndices((maskNumber-1)*d1*d2+1:maskNumber*d1*d2, :) = [baseIndex(:),shiftIndex(:)]; 
    end
end

for k1 = delta1:K1-1
    for k2 = 0:kappa2-1
        maskNumber = maskNumber + 1;
        M_0 = rev2(fft2( M.*circshift( conj(M), [K1-k1 -k2]) ));
        maskComps(:,:, maskNumber) = M_0;
        shiftIndex = circshift(baseIndex, [K1-k1 -k2]);
        sparseIndices((maskNumber-1)*d1*d2+1:maskNumber*d1*d2, :) = [baseIndex(:),shiftIndex(:)]; 
    end
end

for k1 = 0:kappa1-1
    for k2 = delta2:K2-1
        maskNumber = maskNumber + 1;
        M_0 = rev2(fft2( M.*circshift( conj(M), [-k1 K2-k2]) ));
        maskComps(:,:, maskNumber) = M_0;
        shiftIndex = circshift(baseIndex, [-k1 K2-k2]);
        sparseIndices((maskNumber-1)*d1*d2+1:maskNumber*d1*d2, :) = [baseIndex(:),shiftIndex(:)]; 
    end
end

%% 2DWD2 Algorithm

fprintf( ' --- Performing 2DWD2 Algorithm --- \n\n');

tic

% Compute F_D(F_K Y)^T
Y = fft2( permute(fft2(Y), [3 4 1 2]) );

%% Construct vec(X)vec(X)^*

tall_XXstar = zeros(vectLength,1);
maskNumber = 0;

for k1 = 0:kappa1-1
    for k2 = 0:kappa2-1
        maskNumber = maskNumber + 1;
        X_0 = ifft2( Y(:,:,k1+1,k2+1) ./ ( K1*K2*maskComps(:,:,maskNumber) ));
        tall_XXstar((maskNumber-1)*d1*d2+1:maskNumber*d1*d2) = X_0(:);
    end
end

for k1 = delta1:K1-1
    for k2 = delta2:K2-1
        maskNumber = maskNumber + 1;
        X_0 = ifft2( Y(:,:,k1+1,k2+1) ./ ( K1*K2*maskComps(:,:,maskNumber) ));
        tall_XXstar((maskNumber-1)*d1*d2+1:maskNumber*d1*d2) = X_0(:);
    end
end

for k1 = delta1:K1-1
    for k2 = 0:kappa2-1
        maskNumber = maskNumber + 1;
        X_0 = ifft2( Y(:,:,k1+1,k2+1) ./ ( K1*K2*maskComps(:,:,maskNumber) ));
        tall_XXstar((maskNumber-1)*d1*d2+1:maskNumber*d1*d2) = X_0(:);
    end
end

for k1 = 0:kappa1-1
    for k2 = delta2:K2-1
        maskNumber = maskNumber + 1;
        X_0 = ifft2( Y(:,:,k1+1,k2+1) ./ ( K1*K2*maskComps(:,:,maskNumber) ));
        tall_XXstar((maskNumber-1)*d1*d2+1:maskNumber*d1*d2) = X_0(:);
    end
end

XXstar = sparse(sparseIndices(:,1), sparseIndices(:,2), tall_XXstar, d1*d2, d1*d2);

%% Angular Synchronization via Eigenvector Computation

% Hermitian Symmetrize
XXstar = (XXstar + XXstar')/2;

% Recover Magnitude from diag(vec(X)vec(X)^*)
mags = sqrt(diag(XXstar));

% Recover Phase by Computing Eigenvector
diagIndices = find(XXstar);
XXstar(diagIndices) = sign(XXstar(diagIndices));
[X_recon, ~] = eigs(XXstar,1);
X_recon = sign(X_recon);

% Reconstruct Original 2D Signal
X_recon = mags(:) .* X_recon;
X_recon_reshape = reshape(X_recon,d1,d2);

time = toc;

%% Calculate Error in Recovered Signal X

phaseOffset = angle( (X_recon'*X(:))/(X(:)'*X(:)) );
X_phaseCorr = exp( 1i*phaseOffset ) * X_recon;

relError = norm(X(:) - X_phaseCorr,2) / norm(X(:),2);
error_dB = 20 * log10(relError);

%% Print Results

fprintf( ' --------------------------------------------------- \n' );
fprintf( ' Execution Time: %0.3f seconds \n', time );
fprintf( ' --------------------------------------------------- \n\n' );

fprintf( ' Percent Error: %01.2f %% \n', 100*relError );
fprintf( ' Relative Error: %0.3f dB \n', error_dB );

fprintf( '\n --------------------------------------------------- \n' );
fprintf( ' Contact: B. Hutchison   bhutchis@umich.edu \n' );
fprintf( ' --------------------------------------------------- \n\n' );

%% Display Recovered Images

if type == 'b'
    X_phaseCorr = reshape(X_phaseCorr,d1,d2);
    magRecov = abs(full(X_phaseCorr));
    phaseRecov = angle(full(X_phaseCorr));

    figure; subplot(1,2,1);
    imshow(magRecov);
    title('Recovered Mag. Image', 'interpreter', 'latex', 'fontsize', 14);
    xlabel('\textit{Joseph Fourier}', 'interpreter', 'latex', 'fontsize', 14);

    subplot(1,2,2);
    imshow(phaseRecov);
    title('Recovered Phase Image', 'interpreter', 'latex', 'fontsize', 14);
    xlabel('\textit{John Rodenburg}', 'interpreter', 'latex', 'fontsize', 14);
end

%% Functions
function R = rev2(A)
R = circshift( flip(flip(A,1),2), [1 1] );
end
