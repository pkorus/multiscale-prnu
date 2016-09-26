function imsc(x, label, range, structure_field, force_thumbnails)
% imsc(x, [label, range, structure_field, force_thumbnails])
%
% Display an image without axes and print some stats in the title. Images packed in 
% structures and cell arrays will be displayed as grayscale thumbnails.
% 
% Parameters: 
%   - x                - input image / structure / cell array
%   - label            - plot title (%s will be replaced by some basic stats)
%   - range            - range of pixel intensities (for image scaling)
%   - structure_field  - if plotting a cell array containing structures, specify which 
%                        field to plot
%   - force_thumbnails - if plotting a 3-channel image, plot thumbnails of individual
%                        channels instead of a RGB image
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

    if nargin < 2 || isempty(label)
        label = '%s';
    end
    
    if nargin < 3
        range = [];
    end
    
    if strfind(label, '%s')
        format_title = true;
    else
        format_title = false;
    end
    
    if nargin < 5 || isempty(force_thumbnails)
        force_thumbnails = false;
    end
    
    if iscell(x) && numel(x) == 1
        x = x{1};
    end
    
    if ischar(x) && exist(x, 'file')
        x = imread(x);
    end
    
    if isstruct(x)
        fields = fieldnames(x);
        data = cell(1);
        di = 1;
        for i = 1:numel(fields)      
            
            cur_data = x.(fields{i});
            
            if isempty(cur_data)
                continue;
            end
            
            if iscell(cur_data)
                cur_data = cur_data{1};
            end
            
            if isnumeric(cur_data) || islogical(cur_data)
                size_of_data = size(cur_data);
                if size_of_data(1) > 16 && size_of_data(2) > 16
                    data{di} = cur_data;
                    di = di + 1;
                end
            end
        end
        x = data;
    end
    
    if ~iscell(x)
        if size(x, 3) == 1 || (size(x, 3) == 3 && ~force_thumbnails)
            if isempty(range)
                imagesc(x); 
            else
                imagesc(x, range); 
            end
            colors = unique(x);
            if numel(colors) == 2 && size(x, 3) == 1
                if format_title
                    stats = sprintf('avg = %.3f  ccs = %d', mean2(x), max(max(bwlabel(x))));
                    title_nobold(sprintf(label, stats))
                else
                    title_nobold(label);
                end
            else                
                if format_title
                    
                    max_x = max(x(:));
                    if  max_x > 100
                        fstring = '%.0f - %.0f  av = %.0f';
                    elseif max_x > 10
                        fstring = '%.1f - %.1f  av = %.1f';
                    elseif max_x > 1
                        fstring = '%.2f - %.2f  av = %.2f';
                    else
                        fstring = '%.3f - %.3f  av = %.3f';
                    end

                    stats = sprintf(fstring, min(x(:)), max(x(:)), mean2(x));
                    title_nobold(sprintf(label, stats))
                else
                    title_nobold(label);
                end
                
                
            end
        else
            if isempty(range)
                imagesc(generateThumbnails(x, size(x, 3), [], [], true));
            else
                imagesc(generateThumbnails(x, size(x, 3), [], [], false), range);
            end
            if format_title
                stats = sprintf('%d maps along 3-rd dimension', size(x, 3));
                title_nobold(sprintf(label, stats))
            else
                title_nobold(label);
            end
        end
    else
        
        if ~exist('structure_field', 'var')
            thumbs = generateThumbnails(extractCellMaps(x), numel(x), [], [], isempty(range));         
        else
            thumbs = generateThumbnails(extractCellMaps(x, structure_field), numel(x), [], [], isempty(range));
        end
        
        if isempty(range)
            imagesc(thumbs)
        else
            imagesc(thumbs, range)
        end
        
        if format_title
            stats = sprintf('%d maps in the cell array', numel(x));
            title_nobold(sprintf(label, stats))
        else
            title_nobold(label);
        end
    end
    
    cll;    
    
    if size(x,3) == 1
        colormap gray;
    end
end

function title_nobold(string)
    title(string, 'FontWeight', 'Normal', 'FontSize', 10);
end