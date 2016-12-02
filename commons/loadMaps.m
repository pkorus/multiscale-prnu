function response = loadMaps(cam_name, file_name, detector_name, scenario_path, map_type)
% response = loadMaps(cam_name, file_name, detector_name, map_type)
%
% Reads cached response maps for a given test example.
%
% See readme.md for available detector names and map types.

    if ~exist('map_type', 'var')
        map_type = 'tampering_probability';
    end
    
    if ~exist('scenario_path', 'var')
        scenario_path = 'data/tifs-2016-maps';
    end
        
    % Shape of the reliability curve
    w = [30, 2.5];

    data = load(sprintf('%s/%s/%s_%s.mat', scenario_path, cam_name, file_name, detector_name));
    nMaps = numel(data.detection_results.(map_type));
    response = cell(1, nMaps);
    for n = 1:nMaps
        response{n}.candidate = im2double(data.detection_results.(map_type){n});
        response{n}.reliability = 1 - exp(-abs(w(1)*(response{n}.candidate - 0.5).^w(2)));
    end
end