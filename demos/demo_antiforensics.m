
image_name = 'DSC06083'; % 'DSC05748';

gt = imread(sprintf('data/examples/%s_mask.PNG', image_name));
pixmap = imread(sprintf('data/examples/%s.TIF', image_name));
pixmap_org = imread(sprintf('data-input/Sony_A57/pristine/%s.TIF', image_name));

cam = load('data/camera_models/Sony_A57_tiff.mat');

%%

response = detectForgeryPRNUCentral(pixmap, cam.camera_model, 129, ...
    struct('verbose', true, 'stride', 8));

%%

alpha = 0.05;
use_prnu = false;

gtr = repmat(gt, [1 1 3]);

% 
if use_prnu
    pixmap_af = uint8( double(pixmap) .* (1 + alpha *  cam.camera_model.prnu) );
else
    if strcmp(cam.camera_model.stats.denoising, 'bm3d')
        noise = getPRNU(pixmap_org, cam.camera_model.stats.noise_sigma);
    else
        noise = NoiseExtractFromImage(pixmap_org, cam.camera_model.stats.noise_sigma, true);
    end
    noise = ZeroMeanTotal(noise);
    if size(noise,3) > 1
        for s = 1:size(noise,3)
            noise(:,:,s) = WienerInDFT(noise(:,:,s), std2(noise));
        end
    end
    pixmap_af = uint8( double(pixmap) .* (1 + alpha *  noise) );
end

pixmap_af = double(pixmap_af) .* double(gtr) + double(pixmap) .* double(1 - gtr);
pixmap_af = uint8( pixmap_af );

imsc(pixmap_af, sprintf('psnr : %.2f dB', calcPSNR(pixmap, pixmap_af)));

%%

response_af = detectForgeryPRNUCentral(pixmap_af, cam.camera_model, 129, ...
    struct('verbose', true, 'stride', 8));

%%

subplot(2,1,1);
imsc(response, 'forgery')

subplot(2,1,2);
imsc(response_af, 'forgery + anti-forensics')

%%

subplot(2,1,1);
imsc(response.map_cor - response.est_cor, [], [-0.05 0.05]);

subplot(2,1,2);
imsc(response_af.map_cor - response_af.est_cor, [], [-0.05 0.05]);

%% 

subplot(2,1,1);
imsc((response.map_cor - response.est_cor) < - 2 * cam.camera_model.predictors{5}.normfit(2), ...
    'areas with deviation (2*std) from predictions');

subplot(2,1,2);
imsc((response_af.map_cor - response_af.est_cor) < - 2 * cam.camera_model.predictors{5}.normfit(2), ...
    'areas with deviation (2*std) from predictions');