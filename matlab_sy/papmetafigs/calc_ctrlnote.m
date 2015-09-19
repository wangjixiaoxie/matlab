%fxn to calculate mean offset and pvalue for control notes.

function [offsetmn,offsetste, pvl]=calc_ctrlnote(ph,bs)
[plotvls,stats,ctrlvls]=plotinactivfig1v7(ph,bs);

%first do baseline.
indup=find(ctrlvls{1}{1}.drxn==1);
inddn=find(ctrlvls{2}{1}.drxn==2);
for ii=1:2

diff{ii}=[[ctrlvls{ii}{1}.ac(indup)-ctrlvls{ii}{1}.mu(indup)] [ctrlvls{ii}{1}.mu(inddn)-ctrlvls{ii}{1}.ac(inddn)]]

end

combdiff=[diff{1} diff{2}];
offsetmn=nanmean(combdiff);
offsetstd=nanstd(combdiff);
ln=length(find(~isnan(combdiff)));
offsetste=offsetstd./sqrt(ln);

[h,pvl]=ttest(combdiff);