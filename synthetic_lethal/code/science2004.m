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
model = importModel("yeast-GEM.xml");
short_all = SGD(1, 2);
short_all = short_all{1,1};
sys_all = SGD(1, 1);
sys_all = sys_all{1, 1};

doubleKO.genepair1 = {};
doubleKO.genepair2 = {};
doubleKO.growth = {};
doubleKO.exper = {};
tp = 0; fp = 0; fn = 0; tn = 0;
cd ../code
disp('run double gene deletion')
for i = 1:length(data2004)
    g1_sys = short2sys(data2004{i, 1}, short_all, sys_all);
    g2_sys = short2sys(data2004{i, 2}, short_all, sys_all);
    exp = data2004{i, 5};
    [growth, g1, g2] = doubleGD(model, g1_sys, g2_sys);
    if length(growth) == 1
        [tp, tn, fp, fn] = simu_exp(growth, exp, tp, tn, fp, fn);
        doubleKO.genepair1{i} = g1;
        doubleKO.genepair2{i} = g2;
        doubleKO.growth{i} = growth;
        doubleKO.exper{i} = exp;
    end
    disp(i)
end
fprintf('tp=%u\nfp=%u\ntn=%u\nfn=%u\n', tp, fp, tn, fn)
disp('finish')

writecell(doubleKO.genepair1, 'genepair1.xlsx')
writecell(doubleKO.genepair2, 'genepair2.xlsx')
writecell(doubleKO.growth, 'growth.xlsx')
writecell(doubleKO.exper, 'exper.xlsx')
