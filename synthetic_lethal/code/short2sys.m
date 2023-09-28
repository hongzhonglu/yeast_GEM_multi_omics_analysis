%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert gene short name to systematic name
function sys_name = short2sys(short, short_all, sys_all)
    ind = find(strcmp(short, short_all));
    if ~isempty(ind)
        sys_name = sys_all{ind, 1};
    else
        fprintf('%s not in SGD file\n', short)
        sys_name = "nan";
    end
end


