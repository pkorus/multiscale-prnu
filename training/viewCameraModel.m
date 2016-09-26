function viewCameraModel(camera_name, show_scales, print_only, plot_prnu)
% viewCameraModel(camera_name, show_scales, print_only, plot_prnu)
%
% Shows a visualization of multi-scale camera models trained using 'trainCameraModel'.
% 
% Parameters:
%   - camera_name       - name of the camera model (will be sought in ./data/camera_models)
%   - show_scales       - show only plots related to specified scales of analysis, 
%                         e.g., [64, 128] (by default, all scales are included)
%   - print_only        - do not plot figures, just print stats
%   - plot_prnu         - also include a plot of the estimated PRNU (small fragment)
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

    data = load(sprintf('./data/camera_models/%s.mat', camera_name));
    camera_model = data.camera_model;
    
    % Check if the model is complete
    if numel(camera_model.predictors) ~= numel(camera_model.absence_models)
        display('Incomplete model - different lenghts of H1 and H0 models! Training either failed or is incomplete!');
        return;
    end
    
    if nargin < 3 || isempty(print_only)
        print_only = false;
    end
    
    if nargin < 4 || isempty(plot_prnu)
        plot_prnu = false;
    end
    
    % Find matching blocks
    if nargin > 1 && ~isempty(show_scales)
        num_plots = 0;
        for i = 1:numel(camera_model.predictors);
            if any(show_scales == camera_model.absence_models{i}.patch_size)
                num_plots = num_plots + 1;
            end
        end
    else
        num_plots = numel(camera_model.predictors);
    end
    
    if print_only == false
        figure('Name', sprintf('%s - model summary', camera_name)); clf;
        fi = 0;

        glob_maxx = 0;
        for i = 1:numel(camera_model.predictors);
            x = camera_model.predictors{i}.raw_data;
            glob_maxx = max(glob_maxx, max(x));
        end

        for i = 1:numel(camera_model.predictors);

            if nargin > 1 && ~isempty(show_scales) && all(show_scales ~= camera_model.absence_models{i}.patch_size)
                continue;
            end

            fi = fi + 1;

            subplot(4, num_plots, fi);
            data = camera_model.absence_models{i}.raw_data;
            mn = min(data);
            mx = max(data);
            x = linspace(mn, mx,21);
            h = hist(data, x);
            h = h / sum(h) / (x(2)-x(1));    
            handle = bar(x, h, 'hist');
            set(handle, 'FaceColor',[0.85,0.85,1.0]);
            set(handle, 'EdgeColor',[0.85,0.85,1.0]);
            hold on;
            plot(x, normpdf(x, camera_model.absence_models{i}.normfit(1), camera_model.absence_models{i}.normfit(2)), 'r-', 'LineWidth', 1)
            plot([0 0], [0 2*max(h)], 'k--')
            hold off;
            title(sprintf('H0 : %d px window', camera_model.absence_models{i}.patch_size))
            cll;
            xlim([-0.2 0.2])
            ylim([0 1.25*max(h)])

            subplot(4, num_plots, fi + num_plots);
            data = camera_model.predictors{i}.raw_data - camera_model.predictors{i}.raw_predictions;
            mn = min(data);
            mx = max(data);
            x = linspace(mn, mx,65);
            h = hist(data, x);
            h = h / sum(h) / (x(2)-x(1));    
            handle = bar(x, h, 'hist');
            set(handle, 'FaceColor',[0.85,0.85,1.0]);
            set(handle, 'EdgeColor',[0.85,0.85,1.0]);        
            hold on;
            plot(x, normpdf(x, camera_model.predictors{i}.normfit(1), camera_model.predictors{i}.normfit(2)), 'r-', 'LineWidth', 1)
            plot([0 0], [0 2*max(h)], 'k--')
            hold off;
            title(sprintf('H1 : %d px window', camera_model.predictors{i}.patch_size))
            cll;
            xlim([-0.2 0.2])
            ylim([0 1.25*max(h)])

            subplot(4, num_plots, fi + 2*num_plots);
            x = camera_model.predictors{i}.raw_data;
            y = camera_model.predictors{i}.raw_predictions;
            r = rsquare(x, y);
            ranges = cell(1,2);
            ranges{1} = linspace(0,glob_maxx,50);
            ranges{2} = linspace(0,glob_maxx,50);
            n = hist3([x,y], ranges);
            imagesc(-log(n)); cll;
            axis xy;
            hold on;
            plot([0 50], [0 50], 'r:', 'LineWidth', 3);
            hold off;
            title(sprintf('R^2 = %.2f', r))
            xlabel('real correlation')
            ylabel('predicted correlation')
            
            subplot(4, num_plots, fi + 3*num_plots);
            x = camera_model.predictors{i}.raw_data;
            y = camera_model.absence_models{i}.raw_data;

            min_xi = min(min(x), min(y));
            max_xi = max(max(x), max(y));
            xi = linspace(min_xi, max_xi, 35);
            [pi] = hist(x, xi); pi = pi / sum(pi);
            [pj] = hist(y, xi); pj = pj / sum(pj);

            hi = bar(xi, pi, 'hist');
            hold on;
            hj = bar(xi, pj, 'hist');
            hold off;
            set(hi, 'FaceColor', [0.8 1.0, 0.8]);
            set(hj, 'FaceColor', [1.0 0.8, 0.8]);
            cll;
            xlim([min_xi-0.01 max_xi+0.01]);            
            xlabel('real correlation');
            ylabel('frequency');            
        end

        colormap gray;

        if plot_prnu
            figure('Name', sprintf('%s - PRNU signature', camera_name)); clf; colormap gray;
            g_prnu = imcropmiddle(rgb2gray1(camera_model.prnu), [256 256]);
            subplot(1, 2, 1);
            imagesc(g_prnu); title(sprintf('256 x 256 fragment of the PRNU (%s)', strrep(camera_model.name, '_',' ')));
            subplot(1, 2, 2);
            imagesc(abs(fftshift(fft2(g_prnu)))); title('DFT spectrum (magnitude)');
        end
    end
    
    fprintf('\nCamera model:\n');
    fprintf('   name      : %s\n', camera_model.name);
    fprintf('   prnu      : %d x %d x %d\n', size(camera_model.prnu));
    fprintf('   # images  : %d\n', camera_model.stats.prnu_input_count);
    fprintf('   croppig   : %s\n', camera_model.stats.prnu_crop_mode);    
    if isfield(camera_model.stats, 'jpeg_preprocessing')      
        fprintf('   jpeg      : %s\n', camera_model.stats.jpeg_preprocessing)
    else
        fprintf('   jpeg      : %s\n', '-')
    end
    if isfield(camera_model, 'feature_set')      
        fprintf('   features  : %s\n', camera_model.feature_set)
    else
        fprintf('   features  : unspecified (default)\n')
    end
    if isfield(camera_model, 'predictor_type')      
        fprintf('   predictor : %s\n', camera_model.predictor_type)
    else
        fprintf('   predictor : unspecified (linear regression)\n')
    end
    
    fprintf('   denoising : %s with sigma = %.1f\n\n', camera_model.stats.denoising, camera_model.stats.noise_sigma)
    
    fprintf('scale\tR^2\tstd H0\tstd H1 [x 1,000]\n');
    if isfield(camera_model, 'predictors') && isfield(camera_model, 'absence_models')
        for i = 1:numel(camera_model.predictors)
            fprintf('%5d\t', camera_model.predictors{i}.patch_size);
            fprintf('%.2f\t', camera_model.predictors{i}.rsquare);
            fprintf('%5.2f\t', 1000*camera_model.absence_models{i}.normfit(2));
            fprintf('%5.2f', 1000*camera_model.predictors{i}.normfit(2));
            fprintf('\t (predicted from %d patches)', camera_model.predictors{i}.number_of_patches)
            fprintf(' max observed cor = %.2f\n', max(camera_model.predictors{i}.raw_data));
        end
    else
        fprintf('Incomplete model! The H0 or H1 models are missing - please repeat training!\n');
    end    
end