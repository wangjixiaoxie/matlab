function [] = evscreen_auto(batchfile,filetype,pltit,Nmin,thresh);
% usage : evscreen_auto(batchfile,filetype,pltit,Nmin,thresh);
%
% goes through the batch file calculates the spectra
% for each entry, plots it, is considered a song if there are more than
% Nmin segments that
%

if (~exist('pltit'))
    pltit = 0;
end

if (~exist('filetype'))
    filetype = 'ebin0r';
end

Fs = 44100.0;Flw = 300.0;Fhi = 8000.0; %in Hz

if (~exist('Nmin'))
    Nmin = 10;
end

if (~exist('thresh'))
    thresh = 500.0;
end

fkeep = fopen([batchfile,'.keep'],'a');
fdcrd = fopen([batchfile,'.dcrd'],'a');

%figure;m=colormap('gray');m=m(end:-1:1,:);

min_int = 5.0;min_dur = 40.0;smwin=2.0;

fid = fopen(batchfile,'r');
FIRSTWRT=[1,1];
while (1)
    fn  = fgetl(fid);
    if (~ischar(fn))
        break;
    end
    disp(fn);
    if (exist(fn,'file'))
        [rawsong,Fs]=evsoundin('',fn,filetype);
        if (length(rawsong)>0)
            filtsong = bandpass(rawsong,Fs,Flw,Fhi,'hanningfir');
            hlen = floor(Fs*smwin*1.0e-3);
            smooth = conv(ones([hlen,1])./hlen,filtsong.^2);
            [ons,offs]=segment(smooth,Fs,min_int,min_dur,thresh);
            if (length(ons)<Nmin)
                if (FIRSTWRT(1)==1)
                    fprintf(fdcrd,[fn]);
                    FIRSTWRT(1) = 0;
                else
                    fprintf(fdcrd,['\n',fn]);
                end
                if (pltit~=0)
                    hold off;clf;
                    semilogy([1:length(smooth)]/Fs,smooth);
                    hold on;
                    plot([1,length(smooth)]/Fs,[1,1]*thresh,'k--');
                    axis([1/Fs length(smooth)/Fs 1 max(smooth)]);
                    title(fn);
                    drawnow;pause;
                end
            else
                if (FIRSTWRT(2)==1)
                    fprintf(fkeep,[fn]);
                    FIRSTWRT(2) = 0;
                else
                    fprintf(fkeep,['\n',fn]);
                end
            end
        else
            if (FIRSTWRT(1)==1)
                fprintf(fdcrd,[fn]);
                FIRSTWRT(1) = 0;
            else
                fprintf(fdcrd,['\n',fn]);
            end
        end
    end
end
fclose(fid);fclose(fkeep);fclose(fdcrd);
return;
