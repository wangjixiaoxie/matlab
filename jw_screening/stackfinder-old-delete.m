function [scores,spec] = stackfiner(batch,wsize,delay);
%[scores,spec] = stackfinder(batch,wsize,delay);
% find length=delay ms constant-spectrum regions in soundfiles in batch

nfft = wsize;%make obsolete. then should be able to delete
hammy  = hamming(wsize);
ntempl = 2; %would be good to add multiple templates
fs=32000;
Fs=fs; %need this? find/replace

nrep = floor(length(rsong)/wsize)-1;
ntseg = floor(length(rsong)/((delay/1000)*fs))-1; %need this?

scores = zeros([nrep,ntempl]);
spec = zeros([nrep,floor(wsize/2)]);
tspec = zeros([(delay/1000)*fs,size(spec,2)])

filecount=0;

%get filenames from batch
batchid=fopen(batch);
while (~feof(batchid))
	filenm = fscanf(batchid,'%s',1);
	if isempty(filenm)
	        disp('End of soundfiles')
		fclose(batchid);
	        break;
	end

	if exist(filenm)
	        disp(['loading ',filenm]);
	        filecount=filecount+1;
	        [sounddata, fs] = evsoundin('',filenm,'w');
	        sounddata = highpass(sounddata,200,Fs);
	else
	        disp(['file ', filenm,' doesn''t exist on this path?']); continue;
	end
	[s,f,t,p]=spectrogram(sounddata,wsize,.88*wsize,[],32000);
	clear s f t sounddata; pack;
	maxfreq=fix(size(p,1)/2);
	disp(['this file is ',num2str(size(p,2)),' fft segments long.'])


%templates = templates;
% templates are normalized just in case
for jj = 1:ntempl
    if (OLDWAY==0)
        templates(:,jj) = templates(:,jj)-min(templates(:,jj));
        templates(:,jj) = templates(:,jj)./max(templates(:,jj));
    else
        normtmp = templates(:,jj).'*templates(:,jj);
        templates(:,jj) = templates(:,jj)./sqrt(normtmp);
    end
end

for ii = 1:nrep
	ind1 = (ii-1)*nfft+ 1;
	ind2 = ind1 + nfft - 1;
	datchunk = rsong(ind1:ind2) - mean(rsong(ind1:ind2));
	fdatchunk = abs(fft(hammy.*datchunk));

	sp = abs(fdatchunk(1:blen));
    if (USEERROR==0)
        sp(1:6)=0.0;
    else
        %sp(2:end)=sp(1:end-1);
        %sp(1:6)=0.0;
        sp = [zeros([6,1]);sp(6:end-1)];
    end
    if (OLDWAY==1)
        normtmp = sqrt(sp.'*sp);
        sp = sp./normtmp;
    else
        sp = sp-min(sp);
        sp = sp./max(sp);
    end
    for jj = 1:ntempl
        if (OLDWAY==1)
            vals(ii,jj) = acos(sp.'*templates(:,jj));
        else
            vals(ii,jj) = (sp-templates(:,jj)).'*(sp-templates(:,jj));
        end
    end
    spec(ii,:) = sp.';
end
end %for while loop for loading soundfiles from batchfile
return;
