function [] = evscreen(batchfile,filetype,pltspec,specthresh,nfft,OlapFrac,smwin);
%evscreen(batchfile,filetype,pltspec,specthresh,nfft,OlapFrac,smwin);
%
% goes through the batch file calculates the spectra
% for each entry, plots it, can determine it is a song or
% something else
%

if (~exist('specthresh'))
    specthresh = 0.01;
end

if (~exist('filetype'))
    filetype = 'ebin0r';
end

if (~exist('pltspec'))
    pltspec=1;
end

Fs = 44100.0;Flw = 300.0;Fhi = 8000.0; %in Hz
if (~exist('nfft'))
    nfft = 512;
end
spectwin = nfft; 

if (~exist('OlapFrac'))
    noverlap = floor(nfft*0.8);
else
    noverlap = floor(nfft*OlapFrac);
end

if (~exist('smwin'))
    smwin = 2.0; %in ms
end

fkeep = fopen([batchfile,'.keep'],'a');
fdcrd = fopen([batchfile,'.dcrd'],'a');

figure;%m=colormap('gray');m=m(end:-1:1,:);
fid = fopen(batchfile,'r');

while (1)
    fn=fgetl(fid);
    if (~ischar(fn))
        break;
    end
    if (exist(fn,'file'))
        [rawsong,Fs]=evsoundin('',fn,filetype);
    else
        disp(['hey'])
        continue;
    end
    filtsong = bandpass(rawsong,Fs,Flw,Fhi,'hanningfir');
    if (pltspec)
        [sp,f,t] = specgram(filtsong,nfft,Fs,spectwin,noverlap);
        sp = abs(sp);
        p = find(sp<=specthresh);
        sp(p) = specthresh;
    end
    
    hlen = floor(Fs*smwin*1.0e-3);
    smooth = conv(ones([hlen,1])./hlen,filtsong.^2);
    
    if (pltspec)
        ax1=subplot(211);hold off;
        imagesc(t,f,log(abs(sp)));set(gca,'YD','no');
        %colormap(m);colorbar;
        v=axis;axis([v(1:2) 0 9e3]);
        ax2=subplot(212);hold off;linkaxes([ax1,ax2],'x');
    else
        v = [0,length(smooth)/Fs];
    end
    semilogy([1:length(smooth(1:10:end))]*10/Fs,smooth(1:10:end));
    zoom xon;
    vv=axis;axis([v(1:2) 1 vv(4)]);
    title(fn);
    drawnow;
    
    %cmd = input('Keep or toss : [k will keep all else toss] ','s');
    %if (strcmp(cmd,'k'))
    %	fprintf(fkeep,[fn,'\n']);
    %else
    %	fprintf(fdcrd,[fn,'\n']);
    %end
    cmd=menu('What you wanna do?','Keep','Toss','SpectG','Quit');
    if (cmd==4)
        return;
    end
    curfig=gcf;
    while (cmd==3)
        nfig=figure;
        [sp,f,t] = specgram(filtsong,nfft,Fs,spectwin,noverlap);
        sp=abs(sp);p=find(sp<=specthresh);sp(p)=specthresh;
        imagesc(t,f,log(sp));set(gca,'YD','n');ylim([0,1e4]);
        pan xon;zoom xon;
        drawnow;
        cmd=menu('What you wanna do?','Keep','Toss','SpectG','Quit');
        close(nfig);
        if (cmd==4)
            return;
        end
    end
    figure(curfig);
    if (cmd==1)
        fprintf(fkeep,[fn,'\n']);
    elseif (cmd==2)
        fprintf(fdcrd,[fn,'\n']);
    else
        disp(['WHAT?']);
    end
end
fclose(fid);fclose(fkeep);fclose(fdcrd);
return;
