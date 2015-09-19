

%ac710

cd /oriole4/pu34/probein
NT{1}='a'
NT{2}='b'
PRENT{1}=''
PRENT{2}='';
PSNT{1}='';
PSNT{2}='c';
for ii=1:2
[avnac711{ii},t,f]=get_avn('batch10.keep.rand',NT{ii},0.2,0.2,PRENT{ii},PSNT{ii},'obs0'); 
end
% if 'b' is the target note and you only want 'abc' not 'abd' copntext:
%[avn,t,f]=get_avn('batch.train','b',0.2,0.2,'a','c','obs0'); 

cd /oriole4/pu34/lid710
NT{1}='a'
NT{2}='b'
PRENT{1}=''
PRENT{2}='';
PSNT{1}='';
PSNT{2}='c';
for ii=1:2
[avnlid710{ii},t,f]=get_avn('batchlim',NT{ii},0.2,0.2,PRENT{ii},PSNT{ii},'obs0'); 
end

cd /oriole4/pu34/ac714


NT{1}='a'
NT{2}='b'
PRENT{1}=''
PRENT{2}='';
PSNT{1}='';
PSNT{2}='c';
for ii=1:2
[avnac716{ii},t,f]=get_avn('batch16.catch.keep',NT{ii},0.2,0.2,PRENT{ii},PSNT{ii},'obs0'); 
end

cd /oriole4/pu34/lid716
NT{1}='a'
NT{2}='b'
PRENT{1}=''
PRENT{2}='';
PSNT{1}='';
PSNT{2}='c';
for ii=1:2
[avnlid716{ii},t,f]=get_avn('batchlim',NT{ii},0.2,0.2,PRENT{ii},PSNT{ii},'obs0'); 
end
figure
colormap('hot')
figure
ax(1)=subplot(221) 
imagesc(t,f,log(avnac711{1}));syn;ylim([0,1e4]);
mx=max(caxis);
caxis([5 mx])


ax(2)=subplot(222)
imagesc(t,f,log(avnlid710{1}));syn;ylim([0,1e4]);
mx=max(caxis);
caxis([4 mx])

ax(3)=subplot(223)
imagesc(t,f,log(avnac716{1}));syn;ylim([0,1e4]);
mx=max(caxis);
caxis([4 mx])

ax(4)=subplot(224)

imagesc(t,f,log(avnlid716{1}));syn;ylim([0,1e4]);
mx=max(caxis);
caxis([4 mx])

linkaxes(ax);

for ii=1:4
axes(ax(ii));
hold on;
plothorizline([0 6874],[.2 6874],'--','w')
end

for ii=2:4
    axes(ax(ii));
    axis off;
end


indvl=[.075 .075+.032]
for ii=1:2
axes(ax(1));
plotvertline([indvl(ii) 6000],[indvl(ii) 8000],'-','c')
end
