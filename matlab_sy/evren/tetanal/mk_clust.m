function  [denvals,bins,spkampsv]=mk_clust(bt,binwid,TH,BINBND,MAXCNT,...
                                         PLTIT,DOHAD,tet_chans,song_chan);
%[denvals,bins,spkampsv]=mk_clust(bt,binwid,TH,BINBND,MAXCNT,PLTIT,DOHAD,tet_chans,song_chan);
% bt- batch file of cbins
% binwid - width of bins for spk amp histgrams
% TH - spike threshold
% BINBND - boundary of the spkamp density plot bins
% tet_chans - tetrode chans in data file
% DOHAD == 1 means use a hadamard trans on data

if (~exist('DOHAD'))
    DOHAD=0;
end

if (~exist('PLTIT'))
    PLTIT=0;
end

if (~exist('TH'))
	TH=-1000;
else
	if (length(TH)==0)
		TH=-1000;
	end
end
if (~exist('tet_chans'))
	tet_chans=[2:5];
else
	if (length(tet_chans)==0)
		tet_chans=[2:5];
	end
end
if (~exist('song_chan'))
	song_chan=1;
else
	if (length(song_chan)==0)
		song_chan=1;
	end
end

if (length(bt)==0)
	%use the whole directory
	ff=dir('*.cbin');
else
	ff=[];
	fid=fopen(bt,'r');
	cnt=0;
	while (1)
		fn=fgetl(fid);
		if (~ischar(fn))
			break;
		end
		if (exist(fn,'file'))
			cnt=cnt+1;
			ff(cnt).name=fn;
		end
	end
	fclose(fid);
end

BINBND=sort(BINBND);
ll=ceil((BINBND(2)-BINBND(1))./binwid);
comb=[1,2;1,3;1,4;2,3;2,4;3,4];
denvals=[];
if (PLTIT)
    figure;
end
for ii = 1:size(comb,1)
    denvals(ii).denplt=zeros(ll);
    if (PLTIT)
        subplot(3,2,ii);m=colormap('gray');

        cmapl=length(m);
        m=m(end:-1:1,:);
        colormap(m);
        hold on;
        axis([BINBND,BINBND]);
        title([num2str(comb(ii,1)),' vs ',num2str(comb(ii,2))]);
        xlabel(['Chan # ',num2str(comb(ii,1)),' amp']);
        ylabel(['Chan # ',num2str(comb(ii,2)),' amp']);
    end
end
bins=[0:ll-1]*binwid + BINBND(1);

nargout
spkampsv=[];
% cluster all the files in the batch file into 1 plot
for ii = 1:length(ff)
    fn=ff(ii).name;
    disp(fn);
    [data,fs,spkind,spkamp]=tetanal(fn,TH,song_chan,tet_chans);
    if (nargout==3)
        spkampsv=[spkampsv;spkamp];
    end
    if (DOHAD==1)
        spkamp=(hadamard(4)*spkamp.').';
    end
	for jj=1:size(comb,1)
		tmp=denvals(jj).denplt;
		for kk=1:length(spkamp)
			ind1=floor((spkamp(kk,comb(jj,1))-BINBND(1))./binwid)+1;
			ind2=floor((spkamp(kk,comb(jj,2))-BINBND(1))./binwid)+1;
			ind1=min([ll,ind1]);ind1=max([1,ind1]);
			ind2=min([ll,ind2]);ind2=max([1,ind2]);
			tmp(ind2,ind1)=tmp(ind2,ind1)+1;
		end
		denvals(jj).denplt=tmp;
		tmp2=0*tmp;
		pp=find(tmp>0);
		tmp2(pp)=log10(tmp(pp));
        
        if (PLTIT)
            subplot(3,2,jj);imagesc(bins,bins,tmp2);syn;
        end
    end
    if (PLTIT)
        drawnow;
    end
end

if (PLTIT)
    figure;
    for ii = 1:size(comb,1)
        subplot(3,2,ii);

        tmp=denvals(ii).denplt;
        pp2=find(tmp>MAXCNT);
        tmp(pp2)=MAXCNT;
        imagesc(bins,bins,tmp);syn;
        m=colormap('gray');m=m(end:-1:1,:);colormap(m);
        axis([BINBND,BINBND]);
        title([num2str(comb(ii,1)),' vs ',num2str(comb(ii,2))]);
        xlabel(['Chan # ',num2str(comb(ii,1)),' amp']);
        ylabel(['Chan # ',num2str(comb(ii,2)),' amp']);
    end
end
return;
