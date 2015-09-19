%bk41bk14anal.m

%first figure, plot an example song.

figure
ax=gca();

exsong.ax=ax;
exsong.path='/oriole7/dir2/bk67bk44/0927_wnoff'
exsong.fn='bk67bk44_270910_0823.513.cbin';
exsong.bnds=[5.5 6.8]
% 
% exsong.ax=ax;
% exsong.path='/oriole7/dir2/bk67bk44/wnon_913am'
% exsong.fn='bk67bk44_140910_1911.8648.cbin';
% exsong.bnds=[3.2 4.0]


plotcbin(exsong)

% 
% %syntaxanal
% cd /oriole2/bk41bk14/datasum
% load pathvals1-analdata.mat
% for ind=1:length(avls.adjvls{1})
%     
%     for ntvl=1:2
%         ntcnt{ind}{ntvl}=length(avls.adjvls{ntvl}{ind}(:,1))
%     end
% end
% 
% nts=[1 2]
% for ind=1:length(avls.adjvls{1})
%     nt1=ntcnt{ind}{nts(1)}
%     nt2=ntcnt{ind}{nts(2)}
%     sumnts=nt1+nt2
%     ntrat(ind,1)=nt1/sumnts;
%     ntrat(ind,2)=nt2/sumnts;
% end
% avls.ntrat=ntrat;
% save -append pathvals1-analdata.mat avls
%     