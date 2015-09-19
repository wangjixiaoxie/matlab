function []=quickscreenr2(markfile, filetype)
%quick and dirty look at observer data:
%open read each file,and display

if nargin==1
    filetype='obs0r';
end

savefile=[markfile,'.save'];
delfile=[markfile,'.del'];
callfile=[markfile,'.call'];
quietfile=[markfile,'.quiet'];

%open files for saving
if exist(savefile)
    disp([savefile,' already exists']);
    disp(['will append new marks']);
end
save_fid=fopen(savefile,'a');

if exist(delfile)
    disp([delfile,' already exists']);
    disp(['will append new marks']);
end
del_fid=fopen(delfile,'a');

if exist(callfile)
    disp([callfile,' already exists']);
    disp(['will append new marks']);
end
call_fid=fopen(callfile,'a');

if exist(quietfile)
    disp([quietfile,' already exists']);
    disp(['will append new marks']);
end
quiet_fid=fopen(quietfile,'a');

mark='';
counter=0;

%get name of metafile containing file names
meta_fid = -1;
metafile = 0;

disp('select batchfile');
[metafile, pathname]=uigetfile('*','select batchfile')
meta_fid=fopen([pathname, metafile]);
if meta_fid == -1 | metafile == 0
    disp('cannot open file' )
    disp (metafile)
end

while ~strcmp(mark,'q')
counter = counter+1;

%get soundfile name
    soundfile = fscanf(meta_fid,'%s',1);
    if isempty(soundfile)
        disp('End of soundfiles')
        break
    end

%if song exists, get it,
    if exist(soundfile)
        disp('loading file...');
        [sounddata, Fs] = evsoundin('',soundfile,filetype);
    else
        disp(['file ', soundfile,' doesn''t exist on this path?']);
    end

%calculate time indeces for sample points
    time_d=[0:length(sounddata)-1]*1/Fs;  %vector for time axis

%plot file
    if counter==1;
        h_main=figure;
        pos_main=[9 375 560 420];
        fighandle=plot(time_d,sounddata);
        %axeshandle=gca;
        title(soundfile);
        input_flag=1;
    else
        set(fighandle,'YData',sounddata,'XData',time_d);
        title(soundfile);
        input_flag=1;
    end
    while input_flag==1
        figure(h_main);
        set(h_main,'position',pos_main);
        [x,y,mark]=ginput(1);
        mark=char(mark);
        if strcmp(mark,'g')
            %make a spectrogram in a separate window
            h_spect=figure;
            pos_spect=[584 368 560 420];
            set(h_spect,'position',pos_spect);
            [sm,sp,t,f]=evsmooth(sounddata,32000,100);
            imagesc(t,f,log(abs(sp)));set(gca,'YD','n','Ylim',([0 10000]))%'Xlim',xx/1000);
             % specgram(sounddata,512,Fs,[],384);
             mark='';
        elseif strcmp (mark,'p')
            disp([' ',num2str(Fs)]);
            play(sounddata,Fs);
        elseif strcmp(mark,'s')
            fprintf(save_fid,'%s\n',soundfile);input_flag=0;
        elseif strcmp(mark,'d')
            fprintf(del_fid,'%s\n',soundfile);input_flag=0;
        elseif strcmp(mark,'u')
            fprintf(quiet_fid,'%s\n',soundfile);input_flag=0;
        elseif strcmp(mark,'c')
            fprintf(call_fid,'%s\n',soundfile);input_flag=0;
        elseif strcmp(mark,'q')
            disp('quit');input_flag=0;
        end
    end
    if exist('h_spect'); delete(h_spect); clear('h_spect'); end
end
fclose(save_fid);
fclose(del_fid);
fclose(call_fid);
fclose(quiet_fid);

end


