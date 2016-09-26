function colored = colorCodeMap(detection, gt)
% colored = colorCodeMap(detection, gt)
%
    if numel(detection) ~= numel(gt)
        detection = imresize(detection, size(gt), 'nearest');        
    end
    colored = cat(3, gt, detection, detection);
    colored = im2uint8(colored);
end