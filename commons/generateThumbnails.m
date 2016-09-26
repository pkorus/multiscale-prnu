function thumbnails_all = generateThumbnails(images, show_images, margin, thumb_pages, normalize)
% generate_thumbnails(images, show_images, margin, thumb_pages, normalize)
%
% Generates a large image with thumbnails of supplied input images (grayscale only). 
%
% Input params:
% 
%  - images      - 3 dimensional array with input images stacked along the 
%                  3rd dimension
%  - show_images - the number of images per page (scalar);
%                  or the number of images (horizontal/vertical) (2D vec.)
%                  NOTE: if negative, the scanline order changes
%  - margin      - margin around each thumbnail (in pixels)
%  - thumb_pages - limit to the first N pages
%  - normalize   - use 'true' if you want to normalize all images ot [0,1]
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

    if ~exist('margin', 'var') || isempty(margin)
        margin = 1;
    end
    
    if exist('normalize', 'var') && normalize == true
        for i = 1:size(images,3)
            current = images(:,:,i);
            c_min = min(current(:));
            c_max = max(current(:));
            images(:,:,i) = (current - c_min) / (c_max - c_min);
        end
    end

    map_width = size(images,2);
    map_height = size(images,1);

    test_files = size(images,3);

    transpose_mode = false;    

    if ~exist('show_images', 'var'); show_images = min(24, test_files); end;
    
    if any(show_images < 0) 
        transpose_mode = true;
        show_images = abs(show_images);
    end    
    
    if isscalar(show_images)
        images_x = ceil(sqrt(show_images)); 
        images_y = ceil(show_images / images_x); 
        if map_width > map_height && images_x > images_y
            temp = images_x;
            images_x = images_y;
            images_y = temp;
            clear temp
        end
    else
        images_x = show_images(1);
        images_y = show_images(2);
        show_images = prod(show_images);
    end
    if ~exist('thumb_pages', 'var') || isempty(thumb_pages); thumb_pages = ceil(test_files / show_images); end;

    % You can adjust page to "browse" through different pages of results
    last_page = ceil(test_files / show_images);
    thumb_pages = min(last_page, thumb_pages);
    
    wind = [map_height map_width] + margin*2;
    thumbnails_all = zeros(images_y*wind(1), images_x*wind(2), thumb_pages);
    
    for page = 1:thumb_pages
        
        Thumbs = cell(images_y, images_x);
        
        wind = [map_height map_width];
        cnt = 1 + (page-1)*show_images;
        bx = 1;
        by = 1;
        while bx <= images_x && by <= images_y
            if cnt > test_files
                break
            end
            if margin > 0
                Thumbs{by, bx} = padarray(images(:,:,cnt), [margin margin], 1);
            else
                Thumbs{by, bx} = images(:,:,cnt);
            end
            cnt = cnt + 1;
            if transpose_mode == false
                bx = bx + 1;
                if bx > images_x
                    by = by + 1;
                    bx = 1;
                end
            else
                by = by + 1;
                if by > images_y
                    bx = bx + 1;
                    by = 1;
                end
            end
        end
        
        
        wind = wind + margin*2;

        % Collapse cell array into a single image with thumbnails
        thumbnails = ones(size(Thumbs).*wind);
        cnt = 1 + (page-1)*24;
        for bx = 1:images_x
            for by = 1:images_y
                if isempty(Thumbs{by, bx})
                    continue
                end
                thumbnails((by-1)*wind(1)+1:by*wind(1), (bx-1)*wind(2)+1:bx*wind(2)) = Thumbs{by, bx};
                cnt = cnt + 1;
            end
        end
        
        thumbnails_all(:,:,page) = thumbnails;
    end
end