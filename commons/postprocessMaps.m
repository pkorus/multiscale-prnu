function maps = postprocessMaps(maps, fun, field_name)
% maps = postprocessMaps(maps, fun, field_name)
%

    if nargin < 3
        field_name = [];
    end
    
    for i = 1:numel(maps)

        if isempty(field_name)
            for f_name = fieldnames(maps{i})'
                if ~isempty(maps{i}.(f_name{1})) && ~isscalar(maps{i}.(f_name{1})) && ~ischar(maps{i}.(f_name{1}))
                    maps{i}.(f_name{1})= fun(maps{i}.(f_name{1}));
                end
            end
        else
            maps{i}.(field_name)= fun(maps{i}.(field_name));
        end

    end
    
end