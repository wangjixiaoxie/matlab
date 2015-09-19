function [h_ttest,tvl,h_ftest,fvl]=calctftest(vec1,vec2)
[h_ttest,tvl]=ttest2(vec1,vec2);
[h_ftest,fvl]=vartest2(vec1,vec2);