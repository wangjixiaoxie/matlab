function check_notes(batchfile,NOTE,DIVBIN);
% check_notes(batchfile,NOTE,DIVBIN);
%

specthresh = 0.01;
filetype = 'ebin0r';

Fs = 44100.0;Flw = 300.0;Fhi = 8000.0; %in Hz
nfft = 512;
spectwin = nfft; 

noverlap = floor(nfft*0.8);
smwin = 2.0; %in ms

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

    if (exist([fn,'.not.mat'],'file'))
	    load([fn,'.not.mat']);
    else
	    labels=[]; onsets=[]; offsets=[];
    end

    filtsong = bandpass(rawsong,Fs,Flw,Fhi,'hanningfir');
    [sp,f,t] = specgram(filtsong,nfft,Fs,spectwin,noverlap);
    DT = abs(t(2)-t(1));
    sp = abs(sp);p = find(sp<=specthresh);sp(p) = specthresh;
    
    hold off;clf;
    imagesc(t,f,log(abs(sp)));set(gca,'YD','no');
    colorbar;
    v=axis;axis([v(1:2) 0 8e3]);
    title(normstr(fn));
    hold on;grid on;
    plot(v(1:2),[1,1]*f(DIVBIN),'k-');
    
    pp = findstr(labels,NOTE);
    for ii = 1:length(pp)
	    st = floor(onsets(pp(ii))*1e-3/DT)+1;
	    en = floor(offsets(pp(ii))*1e-3/DT)+1;
	    inds = [st:en];
	    tmpsp = sp(1:DIVBIN,inds);
	    [pf,pt]=find(tmpsp==max(max(tmpsp)));

	    tmpsp = sp((DIVBIN+1):end,inds);
	    [pf,pt]=find(tmpsp==max(max(tmpsp)));
	    plot(t(pt+inds(1)-1),f(pf+DIVBIN-1),'ks')
    end

    rdat=readrecf(fn);
    if (length(rdat)>0)
	    tt = rdat.ttimes;
	    for ii = 1:length(tt)
		    plot([1,1]*tt(ii),v(3:4),'k-');
	    end
    end
    zoom xon;
    drawnow;
    
    cmd=menu('What?','Next','Quit');
    if (cmd==2)
	    fclose(fid);
	    return;
    end
end
fclose(fid);
return;
