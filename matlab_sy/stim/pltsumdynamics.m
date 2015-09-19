inds=[17:33]
drxn='up'
figure
amam_stdiff=[]
ampm_stdiff=[]
pmam_stdiff=[]
for ii=inds
    
    amam_stdiff=[amam_stdiff vlout(ii+1).stim(1)-vlout(ii).stim(1)]
    ampm_stdiff=[ampm_stdiff vlout(ii).stim(2)-vlout(ii).stim(1)]
    pmam_stdiff=[pmam_stdiff vlout(ii+1).stim(1)-vlout(ii).stim(2)]
end
if(drxn=='dn')
    amam_stdiff=-amam_stdiff;
    ampm_stdiff=-ampm_stdiff;
    pmam_stdiff=-pmam_stdiff;
end

figure
plot(amam_stdiff, ampm_stdiff,'ko')
hold on;
plot(amam_stdiff, pmam_stdiff,'ko','MarkerFaceColor','k')


