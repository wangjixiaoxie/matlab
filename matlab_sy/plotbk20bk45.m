%raw plot
figure;


path{1}='/oriole4/bk20bk45/probein/'

path{2}='/oriole4/bk20bk45/806_postdata/'

fn{1}='bk20bk45_100708_0711.1064.cbin'
fn{2}='bk20bk45_260808_0805.2570.cbin'

bnds{1}=[2.8335 5.3835]
bnds{2}=[4.73 7.28]

plotlns=1;

axraw(1)=subplot(211)
axraw(2)=subplot(212)

for ii=1:2
    plotcbin(path{ii},fn{ii},bnds{ii},plotlns,axraw(ii));
end

