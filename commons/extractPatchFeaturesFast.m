function features = extractPatchFeaturesFast(patch, std_patch, hpf_std_patch)
%
%  features = extractFeaturesFast(patch, std_patch, hpf_std_patch)
%
%  Extract image features for prediction of PRNU correlation scores.
%   fI  - intensity saturation feature 
%   fT  - texture feature
%   fS  - signal flattening feature
%   fTI - combined texture and intensity feature
%   ... - 2nd order combinations of features
%
% Same as extractFeatures, but requires cached intermediate images for faster computations.
%
% References:
% [1] Chen et al. Determining Image Origin and Integrity using Sensor Noise
%     IEEE Tran. Information Forensics & Security, 2008
%
% -------------------------------------------------------------------------
% This function is a part of multi-scale analysis toolkit available from:
% https://github.com/pkorus/multiscale-prnu
% The code is provided without any warranty or support for educational and 
% research purposes only. See readme.md for more details.
% -------------------------------------------------------------------------
% Written by Paweł Korus, Shenzhen University and AGH University of Science 
%   and Technology
% Version: September 2016
% Contact: pkorus [at] agh [dot] edu [dot] pl
% -------------------------------------------------------------------------

    if size(patch,3) > 1
        patch = rgb2gray(patch);
    end
    
    if ~isfloat(patch)
        patch = double(patch);
    end

    c = 0.03;
    N = numel(patch);
    
    sum2 = @(x) sum(x(:));
    f_att = att(patch, 250, 6);
    f_txt = 1 ./ (1 + hpf_std_patch.^2);
    fI = sum2(f_att) / N;
    fT = sum2(f_txt) / N;     
    fS = mean2(std_patch < c*patch);
    fTI = sum2(f_att .* f_txt) / N;        
    
    features = [1 fI fT fS fTI ...
        fI*fI fI*fT fI*fS fI*fTI fT*fT fT*fS fT*fTI fS*fS fS*fTI fTI*fTI];
end

function y = att(x, crit, tau)
    selector = x > crit;
    y = x / crit;
    y(selector) = exp(-((x(selector) - crit).^2)./tau);
end
