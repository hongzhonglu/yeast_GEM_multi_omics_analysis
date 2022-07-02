function [reaction_expression] = getRxnExp()
fid = fopen('../data/yeast-GEM.txt');
GPR = textscan(fid,'%s %s %s %s %s %s %s','Delimiter','\t','HeaderLines',1);
RxnExp.rxnid = GPR{1};
RxnExp.gpr = GPR{3};
fclose(fid);

for i = 1 : length(RxnExp.gpr)
    temp = split(RxnExp.gpr{i}, ' ');
    % one to one
    if length(temp) == 1
        RxnExp.re{i,1} = RxnExp.gpr{i};
    % isoenzyme
    elseif ismember('or', temp) && ~ismember('and', temp)
        zuhe = strcat('max(', temp{1}, ',', temp{3}, ')');
        if length(temp) > 4
            for j = 5:2:length(temp)
                zuhe = strcat('max(', temp{j}, ',', zuhe, ')');
            end 
        end
        RxnExp.re{i,1} = zuhe;
    % complex
    elseif ismember('and', temp) && ~ismember('or', temp)
        zuhe = strcat('min(', temp{1}, ',', temp{3}, ')');
        if length(temp) > 4
            for j = 5:2:length(temp)
                zuhe = strcat('min(', temp{j}, ',', zuhe, ')');
            end 
        end
        RxnExp.re{i,1} = zuhe;
    % commplex and isoenzyme
    elseif ismember('and', temp) && ismember('or', temp)
        ci = split(RxnExp.gpr{i}, ["(",")"]);
        for c = 2:2:length(ci)
            citemp = split(ci{c}, ' ');
            if length(citemp) > 1
                zuhetemp = strcat('min(', citemp{1}, ',', citemp{3}, ')');
                if length(citemp) > 4
                    for j = 5:2:length(citemp)
                    zuhetemp = strcat('min(', citemp{j}, ',', zuhetemp, ')');
                    end 
                end
                ci{c} = zuhetemp;
            end
        end
        co = 1;
        for f = 1:length(ci)
            if ci{f} ~= "" && ci{f} ~= " or " && length(split(ci{f}, ',')) > 1
                cii{co} = strcat('(',ci{f}, ')');
                co = co + 1;
            elseif ci{f} ~= "" && ci{f} ~= " or " && length(split(ci{f}, ' or ')) > 1
                temp2 = split(ci{f}, ' or ');
                for w = 1:length(temp2)
                    if temp2{w} ~= ""
                        cii{co} = temp2{w};
                        co = co + 1;
                    end
                end
            end
        end
        zuhe = strcat('max(', cii{1}, ',', cii{2}, ')');
        if length(cii) > 2
            for j = 3:length(cii)
                zuhe = strcat('max(', cii{j}, ',', zuhe, ')');
            end 
        end
        RxnExp.re{i,1} = zuhe;
    end
end
reaction_expression = RxnExp.re;

end

