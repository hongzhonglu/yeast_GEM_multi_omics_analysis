%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% predict synthetic lethal
% Data:
% A. H. Tong et al. , Science 294, 2364 (2001)
% A. S. Goehring et al. , Mol. Biol. Cell , 14,1501 (2003)
% K. Kozminski et al. , Mol. Biol. Cell , 14,4958 (2003)
% N. Krogan et al. , Mol. Cell Biol. , 23, 4207 (2003)
% M. Bellaoui et al. , EMBO , 22, 4304 (2003)
% D. Huang et al. , Mol. Cell Biol. , 22, 5076 (2003)
% A. H. Tong et al. , Science (2004)

cd ../data/

disp('prepare data and model')
data2004 = readmatrix("science2004SI.csv", 'OutputType', 'string');
SGDfile = fopen('SGDgeneNames.tsv', 'r');
SGD = textscan(SGDfile, '%s %s', 'Delimiter', '\t');
model = importModel("C:\Users\yuhuzhouye\Desktop\yeast9_w\yesat9_temp\yeast-GEM.xml");
short_all = SGD(1, 2);
short_all = short_all{1,1};
sys_all = SGD(1, 1);
sys_all = sys_all{1, 1};


doubleKO.genepair1 = [];
doubleKO.genepair2 = [];
doubleKO.growth = [];
doubleKO.exper = [];

disp('run double gene deletion')
for i = 1:length(data2004)
    g1_sys = short2sys(data2004{i, 1}, short_all, sys_all);
    g2_sys = short2sys(data2004{i, 2}, short_all, sys_all);
    exp = data2004{i, 5};
    [growth, g1, g2] = doubleGD(model, g1_sys, g2_sys);
    doubleKO.genepair1(i, 1) = g1;
    doubleKO.genepair2{i, 1} = g2;
    doubleKO.growth{i, 1} = growth;
    doubleKO.exper{i, 1} = exp;
    disp(i)
end
disp('finish')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert gene short name to systematic name
function sys_name = short2sys(short, short_all, sys_all)
    ind = find(strcmp(short, short_all));
    if ~isempty(ind)
        sys_name = sys_all{ind, 1};
    else
        fprintf('%s not in SGD file\n', short)
        sys_name = 'nan';
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run double gene deletion
% gene is in systematic name
function [growth, g1, g2] = doubleGD(model, gene1, gene2)
    ind1 = find(strcmp(gene1, model.genes));
    ind2 = find(strcmp(gene2, model.genes));
    if ~isempty(ind1) && ~isempty(ind2)
        fprintf('deleting %s and %s\n', gene1, gene2)
        [~, grRateKO, ~] = doubleGeneDeletion(model,'FBA', model.genes(ind1), model.genes(ind2));
        growth = grRateKO;
        g1 = gene1;
        g2 = gene2;
        disp(growth)
    else
        fprintf('%s or %s not in model\n', gene1, gene2)
        growth = 'nan';
        g1 = 'nan';
        g2 = 'nan';
    end
end
