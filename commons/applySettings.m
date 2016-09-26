function settings = applySettings(defaults, user, supported_fields)
% settings = applySettings(defaults, user, supported_fields)

    settings = struct();
    names = fieldnames(defaults);
    
    if nargin < 3
        supported_fields = [];
    end
    
    if ~isempty(supported_fields)
        set_fields = fieldnames(user);
        unsup_fields = setdiff(set_fields, supported_fields);
        if numel(unsup_fields) > 0
            throw(MException('prnu:UnknownParameter', sprintf('Unsupported parameters: %s', strjoin(unsup_fields))))
        end
    end
    
    for i = 1:numel(names)
        if isfield(user, names{i})
            settings.(names{i}) = user.(names{i});
        else
            settings.(names{i}) = defaults.(names{i});
        end                
    end

end