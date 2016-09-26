function score = scoreLocalization(b_map, gt)
%
% score = scoreLocalization(b_map, gt)
%
% Computes popular tampering localization quality metrics.
% 
% The output structure contains the following fields:
%
%         tpr: true postive rate
%         tnr: true netagive rate
%    accuracy: balanced accuracy (tpr+tnr)/2
%   precision: precision
%      recall: recall
%          f1: harmonic mean of precision and recall
%
    
    if size(b_map,1) ~= size(gt,1) || size(b_map,2) ~= size(gt,2)
        b_map = imresize(b_map, size(gt), 'nearest');
    end

    score = struct();
    sum2 = @(x) sum(x(:));
    
    TP = sum2(b_map & gt);          % detected forged pixels
    TN = sum2((1-b_map) & (1-gt));  % correctly detected genuine pixels

    precision = TP / sum2(b_map);
    if isnan(precision); precision = 0; end;
    recall = TP / sum2(gt);

    f1 = 2 * (precision * recall) / (precision + recall);
    if isnan(f1); f1 = 0; end;
    tpr = recall;
    tnr = TN/sum2(1-gt);
    accuracy = (tpr + tnr)/2;

    score.tpr = single(tpr);
    score.tnr = single(tnr);
    score.accuracy = single(accuracy);
    score.precision = single(precision);
    score.recall = single(recall);
    score.f1 = single(f1);

end