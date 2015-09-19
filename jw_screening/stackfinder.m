function [scores,spec] = stackfiner(batch,wsize,delay);
%[scores,spec] = stackfinder(batch,wsize,delay);
% find length=delay ms constant-spectrum regions in soundfiles in batch

fs=32000;
ntempl = 2; %would be good to add multiple templates
Fs=fs; %need this? find/replace
overlap=.8;
delay=round(delay*(fs/1000)/((1-overlap)*wsize))


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
	        sounddata = highpass(sounddata,200,fs);
		disp('filtered.')
	else
	        disp(['file ', filenm,' doesn''t exist on this path?']); continue;
	end
	[s,f,t,p]=spectrogram(sounddata,wsize,fix(overlap*wsize),2*floor(wsize/2),fs);
	clear s f sounddata;
	nrep=size(p,2);
	maxfreq=fix(size(p,1)/2);
	disp(['this file is ',num2str(nrep),' fft segments long.'])
	spec=p;
	scores = zeros([nrep,ntempl+1]);

for ii = 1:nrep

        %sp = [zeros([6,1]);sp(6:end-1)]; %skipping this due to use of highpass above
        skip=0;
	if delay > ii-1
		skip=1; %change to disabling template 1
	end
	if delay > (nrep-ii)-1
		if ii==nrep
			save([filenm,'.stackscore'],'scores');	
		end
		skip=2; %change to disabling template 2
	end
	sp=p(:,ii);
	sp = sp-min(sp);
        sp = sp./max(sp);
%	templates(:,1) = p(:,ii-delay);
%	templates(:,2) = p(:,ii+delay);

	for jj = 1:ntempl
		if skip==1;if jj==1;continue;end;else templates(:,1) = p(:,ii-delay);end;
			if skip==2;if jj==2;continue;end;else templates(:,2) = p(:,ii+delay);end
%		templates(:,1) = p(:,ii-delay);
%	        templates(:,2) = p(:,ii+delay);
		templates(:,jj) = templates(:,jj)-min(templates(:,jj));
		templates(:,jj) = templates(:,jj)./max(templates(:,jj));
           	scores(ii,jj) = (sp-templates(:,jj)).'*(sp-templates(:,jj));
	end
    
	if ii==nrep %currently useless but when template disableing is added it will run. In which case take out the duplicate if loop above
		save('test','scores');
	end
end %for ii loop
end %for while loop for loading soundfiles from batchfile
scores(:,3)=t(1,:);
return;
