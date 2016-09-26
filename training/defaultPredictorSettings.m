function o = defaultPredictorSettings()
% o = defaultPredictorSettings()
%
%     patches_per_image: + integer, default 500
%                  seed: + integer, default 1234
%           noise_sigma: float, default 3
%          noise_method: 'wavelet' or 'bm3d'
%          jpeg_quality: number of 2-D vector (range), e.g., [75-100]
%           feature_set: 0 - standard 15-D feature set
%                        1 - standard + extra JPEG feature
%        predictor_type: 'ls' - least square
%                        'nn' - neural network
%                nIters: +integer, default 1

    o = struct();
    o.patches_per_image = 500;
    o.seed = 1234;
    o.noise_sigma = 3;
    o.noise_method = 'wavelet';
    o.jpeg_quality = [];
    o.feature_set = 0;
    o.predictor_type = 'lr';
    o.nIters = 1;    

end