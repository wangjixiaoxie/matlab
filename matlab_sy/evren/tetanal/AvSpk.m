%get the average spike wavform for the cluster given in
%CURRCLUSTERBOX

if (~exist('clustdat'))
	return;
end

curclust = str2num(get(CURRCLUSTERBOX,'String'));
if ((curclust<1)|(curclust>length(clustdat)))
	return;
end

if (~exist('AVSPKPLT'))
	AVSPKPLT=figure;
	set(AVSPKPLT,'Tag','AverageSpkPlt','Position',[711,41,560,420],...
	                      'DeleteFcn','AVSPKPLTDELETE');

	ax=zeros([4,1]);
	for ii = 1:4
		ax(ii)=subplot(2,2,ii);
	end
	linkaxes(ax);
else
	figure(AVSPKPLT);
	for ii=1:4
		subplot(2,2,ii);hold on;
	end
end

clr    = clustdat(curclust).clr;
s_inds = spkind(clustdat(curclust).spikes);

% how far back and forward to go 
av_rng  = [-5,10]*1e-3;
av_inds = [floor(av_rng(1)*fs):ceil(av_rng(2)*fs)];
avspk   = zeros([length(av_inds),4]);
nn=0;
for ii=1:length(s_inds)
	%avoid the few spikes within av_rng(1) of begining and
	%av_rng(2) of end of data file
	if ((s_inds(ii)+av_inds(1)>0)&((s_inds(ii)+av_inds(end))<size(data,1)))
		nn=nn+1;
		for jj = 1:4
			avspk(:,jj) = avspk(:,jj) + ...
			         data(s_inds(ii)+av_inds,tet_chans(jj));
		end
	end
end
avspk=avspk/nn;
for ii=1:4
	subplot(2,2,ii);
	plot(av_inds*1e3/fs,avspk(:,ii),[clr,'.-']);
	title(['Chan ',num2str(ii)]);
    tmplim  = ylim;
    tmplim2 = [min(min(avspk))*1.1,max(max(avspk))*1.1];
    ylim([min([tmplim(1),tmplim2(1)]),max([tmplim(2),tmplim2(2)])]);
	grid on;
end
zoom on;
