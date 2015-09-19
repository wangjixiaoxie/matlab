function[threshval,threshind]=calc_escapethresh(vals,ps)
threshind=find(vals(:,3)==1);
threshval=median(vals(threshind,2));

