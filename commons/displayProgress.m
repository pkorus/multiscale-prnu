function displayProgress(label, progress, count)

    if progress > 0
	for pp = 1:(5+numel(label)+count+5);
	    fprintf('\b');
	end
    else
	fprintf('\n');
    end
    fprintf('%s : [', label);
    for pp = 1:progress
	fprintf('#');
    end
    for pp = progress+1:count
	fprintf(' ');
    end
    fprintf('] %3d%%', round(progress/count*100)); 
            
end