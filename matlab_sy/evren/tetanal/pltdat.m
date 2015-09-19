if (~exist('DATPLT'))
	DATPLT=figure;
	set(gcf,'Pos',[5,40,664,909]);
	set(gcf,'Tag','RawDataFig');
	set(gcf,'DeleteFcn','PLTDATDELETE');
	MOVELEFTBTN=uicontrol('Style','pushbutton','Position',[90,40,50,20],...
	    'String','<-','CallBack','MoveLeft');
	MOVERGHTBTN=uicontrol('Style','pushbutton','Position',[150,40,50,20],...
	    'String','->','CallBack','MoveRght');
    PLTSCATBTN=uicontrol('Style','pushbutton','Position',[210,40,75,20],...
	    'String','Scatter','CallBack','pltscat');
else
	figure(DATPLT);
	for ii = 1:5
		subplot(5,1,ii);
		cla;
	end
end
Nusechans=length(tet_chans);
ax=zeros([Nusechans+1,1]);
ax(1)=subplot(Nusechans+1,1,1);
[sm,sp,t,f]=evsmooth(data(:,song_chan),fs,10);
imagesc(t,f,log(abs(sp)));syn;ylim([0,1e4]);
rd=readrecf(fn);
[pth,nm,ext]=fileparts(rd.outfile);
title([normstr(fn),' - ',rd.outfile],'Interpreter','none');

for ii = 1:Nusechans
	ax(ii+1)=subplot(Nusechans+1,1,ii+1);
	plot([1:length(data)]/fs,data(:,tet_chans(ii)));
	grid on;
end

%[trig,fs]=evsoundin('',fn,'obs5');
%if (DOTRIGPLT)
%    hold on;TRIGPLT=plot([1:length(trig)]/fs,trig*1e4/max(trig),'k');
%end

linkaxes(ax,'x');
pan xon;zoom xon;

set(DATPLT,'Toolbar','figure');
