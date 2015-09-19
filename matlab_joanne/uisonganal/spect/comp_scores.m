function [ratscore1,ratconf]=comp_scores(scores, ind1, ind2,  v1, v2);



diff_flag = 0; %if 1 use differences instead of ratios

aut1=v1(ind1);
aut2=v2(ind2);

aut_ratio1=aut1./aut2;
aut_ratio2=aut2./aut1;
aut_diff=aut1-aut2;

ratscore1=2*ones(1,length(aut_ratio1));
ratscore1_ind=find(aut_ratio1>=1);
ratscore1(ratscore1_ind)=ones(1,length(ratscore1_ind));
ratconf=max(aut_ratio1,aut_ratio2);


if diff_flag==1
  diffscore1=2*ones(1,length(aut_ratio1));
  diffscore1_ind=find(aut_diff>=0);
  diffscore1(diffscore1_ind)=ones(1,length(diffscore1_ind));
  diffconf=abs(aut_diff);
  ratscore1=diffscore1;
  ratconf=diffconf;
end

%figure
%subplot(2,1,1);
%plot(ratconf,ratscore1,'+');
%subplot(2,1,2);
%plot(diffconf,diffscore1,'+');


%calcualte and display concordance between blind scores and auto scores
figure
subplot(2,1,1)
agreement=abs(scores(1,:)-ratscore1);
agreement=abs(1-agreement);
[x,cidx]=sort(ratconf);
cidx=fliplr(cidx);
ragreement=agreement(cidx);
plot(100*ragreement,'+');

flen=20
filtper=100/flen*ones(1,flen);
ragreement=conv(filtper,ragreement);
keep_idx=[flen:length(ragreement)-flen];
ragreement=ragreement(keep_idx);
rflen=(flen/2);
plot_idx=keep_idx-rflen;
hold on
plot(plot_idx,ragreement,'m');
plot([0,length(ind1)],[50,50],'r:');
ylabel('% agreement between auto and blind scores');
xlabel('confidence (rank order of auto measure)');



%plot performance of blind and auto measures

subplot(2,1,2);

%auto measures
agreement=abs(ratscore1-2);
[x,cidx]=sort(ratconf);
cidx=fliplr(cidx);
ragreement=agreement(cidx);
rscores=scores(1,:);
rscores=rscores(cidx);
rscores=abs(rscores-2);
plot(100*ragreement,'+');
hold on
plot(100*rscores,'mo');

flen=20
filtper=100/flen*ones(1,flen);
ragreement=conv(filtper,ragreement);
keep_idx=[flen:length(ragreement)-flen];
ragreement=ragreement(keep_idx);
rflen=(flen/2);
plot_idx=keep_idx-rflen;
hold on
plot(plot_idx,ragreement);
plot([0,length(ind1)],[50,50],'r:');
ylabel('% of first arg ranked as less noisy');
xlabel('confidence (rank order of auto measure)');

%blind measures ploted as function of auto confidence

agreement=abs(scores(1,:)-2);
[x,cidx]=sort(ratconf);
cidx=fliplr(cidx);
ragreement=agreement(cidx);

flen=20
filtper=100/flen*ones(1,flen);
ragreement=conv(filtper,ragreement);
keep_idx=[flen:length(ragreement)-flen];
ragreement=ragreement(keep_idx);
rflen=(flen/2);
plot_idx=keep_idx-rflen;
hold on
plot(plot_idx,ragreement,':');


%blind measures ploted as function of blind confidence

agreement=abs(scores(1,:)-2);
[x,cidx]=sort(scores(2,:));
cidx=(cidx);
ragreement=agreement(cidx);
hold on
flen=20
filtper=100/flen*ones(1,flen);
ragreement=conv(filtper,ragreement);
keep_idx=[flen:length(ragreement)-flen];
ragreement=ragreement(keep_idx);
rflen=(flen/2);
plot_idx=keep_idx-rflen;
hold on
plot(plot_idx,ragreement,'--g');


title('- aut vs aut_c; .. blind vs aut_c; -- blind vs blind_c');

subplot(2,1,1)







