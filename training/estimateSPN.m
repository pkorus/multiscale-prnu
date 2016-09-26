function [RP, LP, uImages] = estimateSPN(Images, spn_type, denoiser, sigma, sepChannels) 
% [RP, LP, uImages] = estimateSPN(Images, spn_type, denoiser, sigma, sepChannels) 
%
% Estimates PRNU pattern from given images. This function is a
% generalization of 'getFingerprint' from [1]. New features include:
%
%  - support for both MLE estimation and phase SPN estimation [2]
%  - support for both mihcak and bm3d denoising
%
% In case of MLE estimation, the PRNU will be post-processed by the
% standard 'zero-mean' procedure and by Wiener filtering in DFT domain.
%
% Parameters:
%
%  - Images       - list of file names
%  - spn_type     - 'mle' or 'phase'
%  - denoiser     - 'mihcak' or 'bm3d'
%  - sigma        - noise standard deviation
%  - sepChannels  - whether to keep separate RGB channels or not
%
% The recommended way of building camera models is through the 
% trainCameraModel() function.
%
% References:
%
% [1] http://dde.binghamton.edu/download/camera_fingerprint/
%
% [2] X. Kang et al., Enhancing Source Camera Identification Performance 
%     with a Camera Reference Phase Sensor Pattern Noise, IEEE Tran. on
%     Information Forensics and Security, doi: 10.1109/TIFS.2011.2168214
% -------------------------------------------------------------------------
% This function is a part of multi-scale analysis toolkit available from:
% https://github.com/pkorus/multiscale-prnu-localization-toolbox
% The code is provided without any warranty or support for educational and 
% research purposes only. See readme.md for more details.
% -------------------------------------------------------------------------
% Written by Paweł Korus, Shenzhen University and AGH University of Science 
%   and Technology
% Version: September 2016
% Contact: pkorus [at] agh [dot] edu [dot] pl
% -------------------------------------------------------------------------

    nImages = numel(Images);
    
    if nImages == 0
        throw(MException('prnu:InvalidParameter', 'Empty image list!'));
    end
    
    if ~exist('sigma', 'var')
        sigma = 3.0;
    end
    
    if ~exist('sepChannels', 'var')
        sepChannels = true;
    end
    
    fprintf('PRNU estimation: %s + %s(%.1f) denoiser\n\n', spn_type, denoiser, sigma);
    
    last_progress = -1;
    uImages = 0;
    for n = 1:nImages       
        progress = round(100*(n/nImages)/2);
        
        % Display progress bar
        if progress ~= last_progress
            displayProgress(sprintf('Extracting noise residuals : %4d / %4d', n, nImages), progress, 50);
        end
        last_progress = progress;  
        
        if isstruct(Images)
            X = imread(Images(n).name);
        else
            X = imread(Images{n});
        end
        X = double255(X);
        
        if n == 1
            
            [nRows,nCols,nChannels] = size(X);
            
            if nChannels ~= 3 
                continue;
            end
            
            %%%  Initialize sums
            RPsum = cell(1, 3);
            NN = cell(1, 3);
            
            for j = 1:3
                RPsum{j} = zeros(nRows, nCols, 'single');   
                NN{j} = zeros(nRows, nCols, 'single');        	% number of additions to each pixel for RPsum
            end
            
        else
            s = size(X);
            if length(size(X))~=3, 
                continue;                           % only color images will be used 
            end
            if any([nRows,nCols,nChannels] ~= s)
                continue;                           % only same size images will be used 
            end
        end
        
        if strcmpi(spn_type, 'mle')
            for j = 1:3
                ImNoise = single(extractNoise(X(:,:,j), denoiser, sigma)); 
                Inten = single(IntenScale(X(:,:,j))) .* Saturation(X(:,:,j));    % zeros for saturated pixels
                RPsum{j} = RPsum{j} + ImNoise .* Inten;   	% weighted average of ImNoise (weighted by Inten)
                NN{j} = NN{j} + Inten.^2;
            end
        elseif strcmpi(spn_type, 'phase')
            for j = 1:3
                ImNoise = single(extractNoise(X(:,:,j), denoiser, sigma));                 
                W = fft2(ImNoise);
                W = W ./ abs(W);                
                RPsum{j} = RPsum{j} + W;
            end
        else
            throw(MException('prnu:InvalidParameter', 'Unsupported SPN type! Use either: mle or phase'));            
        end
        
        uImages = uImages + 1;        
    end
        
    if uImages == 0
        throw(MException('prnu:InvalidImages', 'None of the images was a color image in landscape orientation!'));
    else
        fprintf('Successfully estimated PRNU from %d images\n', uImages);
    end
    
    RP = zeros(nRows, nCols, nChannels);
    if strcmpi(spn_type, 'mle')
        for j = 1:3
            RP(:,:,j) = RPsum{j}./(NN{j}+1);
        end
    elseif strcmpi(spn_type, 'phase')
        for j = 1:3
            RP(:,:,j) = real(ifft2(RPsum{j}./uImages));
        end
    end    
    
    % Remove linear pattern and keep its parameters
    [RP, LP] = ZeroMeanTotal(RP);
    
    % Process color channels according to the requirements
    if sepChannels
        if strcmpi(spn_type, 'mle')
            for i = 1:3
                RP(:,:,i) = WienerInDFT(RP(:,:,i),std2(RP(:,:,i)));
            end
        end
    else
        RP = rgb2gray1(RP);
        if strcmpi(spn_type, 'mle')
            RP = WienerInDFT(RP, std2(RP));
        end
    end
    
    % Reduce double to single precision
    RP = single(RP);        
end    

function noise = extractNoise(image, denoiser, sigma)

    if strcmpi(denoiser, 'bm3d')
        [~, noise] = BM3D([], image, sigma);
        noise = image - 255*noise;
    elseif strcmpi(denoiser, 'mihcak')
        L = 4;
        qmf = MakeONFilter('Daubechies',8);
        noise = NoiseExtract(image, qmf, sigma, L);
    else
        throw(MException('prnu:InvalidParameter', 'Unsupported denoiser! Use either: mihcak or bm3d'));
    end
end

function X = double255(X)
% convert to double ranging from 0 to 255
    datatype = class(X);
    switch datatype,                % convert to [0,255]
        case 'uint8',  X = double(X);
        case 'uint16', X = double(X)/65535*255;  
    end
end