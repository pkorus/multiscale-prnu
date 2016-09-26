function y = extractField(c, f)
    y = zeros(size(c));
    for i = 1:numel(c)
        y(i) = c{i}.(f);
    end
end