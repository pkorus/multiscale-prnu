function thresholds = getThresholds(num)
    thresholds = linspace(0,1,num+2); 
    thresholds = thresholds(2:end-1);
end