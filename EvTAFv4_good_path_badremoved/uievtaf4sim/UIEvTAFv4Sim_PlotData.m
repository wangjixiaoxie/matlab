function UIEvTAFv4Sim_PlotData(hObject,handles,REPLOTSPECGRAM);
%uievtafsim_plotdata(hObject,handles);

set(handles.FileNumBox,'String',[num2str(handles.Nfile),'/',num2str(length(handles.files))]);
t=handles.t;fs=handles.fs;dat=handles.dat;sp=handles.sp;labels=handles.labels;
onsets=handles.onsets;offsets=handles.offsets;

fn=handles.files{handles.Nfile};

if (REPLOTSPECGRAM==1)
    
    [dat,fs]=evsoundin('',fn,handles.CHANSPEC);

    [sm,sp,t,f]=evsmooth(dat,fs,0);
    sp=abs(sp);
    if (~exist([fn,'.not.mat'],'file'))
        min_int = 5;
        min_dur = 30;
        sm_win=2;
        [onsets,offsets]=SegmentNotes(sm,fs,min_int,min_dur,handles.SEGTH);
        onsets = onsets*1e3;offsets=offsets*1e3;
        labels = char(ones([1,length(onsets)])*fix('-'));
    else
        load([fn,'.not.mat']);
        %if (handles.SEGTH~=threshold)
        %    [onsets,offsets]=SegmentNotes(sm,Fs,min_int,min_dur,threshold);
        %    labels = char(ones([1,length(onsets)])*48);
        %end
    end

    axes(handles.SpecGramAxes);hold off;
    imagesc(t,f,log(abs(sp)));
    set(gca,'YD','n');
    set(gca,'YTick',[0:2e3:1e4]);
    set(gca,'YTickLabel',{'0','2','4','6','8','10'});
    title(fn,'Interpreter','none');xlim([0,t(end)]);ylim([0,1e4]);
end

handles.labels=labels;
handles.SPMax=max(max(log(sp)));
handles.SPMin=min(min(log(sp)));
temp1=get(handles.SPMinLevel,'Value');temp2=get(handles.SPMaxLevel,'Value');
mn=handles.SPMin;mx=handles.SPMax;
caxis([temp1*(mx-mn)+mn,temp2*mx]);

handles.t=t;
handles.fs=fs;
handles.dat=dat;
guidata(hObject,handles);

PlotTafValsv4(hObject,handles,REPLOTSPECGRAM);
handles=guidata(hObject);

axes(handles.LabelAxes);cla;hold off;
guidata(hObject,handles);
if (exist('handles.ActTrigTimes'))
    delete(handles.ActTrigTimes);
    handles.ActTrigTimes=[];
end
set(handles.LabelAxes,'XTick',[],'YTick',[]);

if (length(labels)>0)
    handles.LabelHandl=text((onsets+offsets).'*5e-4,0*onsets.',labels.');hold on;
else
    handles.LabelHandl=[];
end


rdata=readrecf(fn);
if (exist('rdata'))
    if (length(rdata.ttimes)>0)
        handles.ActTrigTimes=plot(rdata.ttimes*1e-3,0*rdata.ttimes,'rv');hold on;
    end
end

trigtimes = handles.TrigT;
tmp_plt=plot(trigtimes,0*trigtimes-0.5,'b^');
axis([0,t(end),-1,1]);
handles.TRIGPLT=tmp_plt;

linkaxes([handles.SpecGramAxes,handles.ValAxes,handles.LabelAxes],'x');

handles.tlim=[0,t(end)];
handles.sp=sp;
handles.fs=fs;
handles.onsets=onsets;
handles.offsets=offsets;
handles.TrigT=trigtimes;
guidata(hObject,handles);

return;
