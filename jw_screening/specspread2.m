function [specmat,matdev,avgspec]=specspread2(batch)
%[specimg]=specspread1(batch)
%in specspread2, information across segments is being discarded for speed. stats only poss using specspread1.
%for specspread2 to work, had to transpose columns of p during the addition.. no change in final output. why?

%using lag=0. assuming Fs=44100 and w type file for evsoundin.
%for ZF y75bu48... 100ms and 2-300ms may be interesting lags, based on 21d post-deaf. 20-21seg/100ms.
lag=0;
Fs=44100;
prevz=0; %counter for concatinating the z-planes in specmat.
filecount=0;

%%from single-file version. hard to init with multiple files.
%init matrix here. spectrogram outputs 513 f bins for 512 window. NOT TRUE??
%maxfreq=fix(size(p,1)/2);
%specmat((maxfreq)+1,:,:)=0; %the first row will be the unscaled spectrum of that segment
%specmat(:,maxfreq+1,:)=0; %the first column will be the scalar (saving explicitly allows for non-zero lag)
%specmat(:,:,size(p,2))=0; %without knowing how many files and segs/file, how to init the z-dimension?

%get filenames from batch
batchid=fopen(batch);
%%%whileloop.. indent..%%%%
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
	randy=fix(rand*20);
	[sounddata, Fs] = evsoundin('',filenm,'w');
	sounddata = highpass(sounddata,200,Fs);

else
	disp(['file ', filenm,' doesn''t exist on this path?']); continue;
end

[s,f,t,p]=spectrogram(sounddata,512,300,[],44100); %sounddata,512,300,[],44100 means 212samples per segment.
maxfreq=fix(size(p,1)/2);
disp(['this file is ',num2str(size(p,2)),' fft segments long.'])
clear s f t sounddata;
if filecount==1
	specmat(maxfreq+1,maxfreq+1)=0;
	matdev(maxfreq+1,maxfreq+1)=0;
end

for z=1:size(p,2)
	%%progress report:
	if fix(z/500)==z/500	
		disp(['       ',num2str(z),'  ',num2str(prevz),'  ',num2str(size(specmat,2)),'  ',num2str(size(specmat,1))])
	end
	
	specmat(1,2:end)=specmat(1,2:end)+p(1:maxfreq,z)'; 
	specmat(2:end,1)=p(1:maxfreq,z);
	
	for i=2:size(specmat,1)
		if z+lag > size(p,2) | z+lag <= 0
	%		specmat(2:end,2:end)=0;
			break;
		end
		toadd=(p(1:maxfreq,z+lag)'*specmat(i,1));
		if z+prevz==1
			avgest=toadd; %this is to avoid /0 err. for first window, we can skip the addition of deviation totals.
		else
			avgest=((specmat(i,2:end)./(prevz+z-1)));
		end
		% standard deviation estimation:
		if filecount ~=1 %to exclude first few windows, better std estimate. and can get rid of if-loop above?
		matdev(i,2:end)=matdev(i,2:end)+(toadd-avgest).^2; %in specspread1, matdev was almost always about .88 to .98 times the real std when not excluding first windows. Oops, had the "randy"-based sampler running too.
		else
		varz=z;
		end
		%% end var data collection
		specmat(i,2:end)=specmat(i,2:end)+toadd;
	       	
	end
end
%% important. here, create/concat a txt file saying what file gets what "z" range in specmat.
% or create var here and save at the end.
filehist(filecount,:)={filenm,'starts after',num2str(prevz),'goes till',num2str(prevz+z)};
prevz=prevz+z;

end %%%end whileloop.. indent..%%%%

specmat=(specmat./prevz);
avgspec=specmat(1,2:end);
specmat=specmat(2:end,2:end);
%specimg=flipud(specimg(2:end,:));
matdev=sqrt(matdev(2:end,2:end)./(prevz-varz-1)); %to use matdev to report stdev, comment out next line. Checked with highflat.. this is a good std estimate.

%matdev=matdev./specmat; %this line to use matdev to report CV.

%specmat(5000:100:7000)
%matdev(5000:100:7000)
%jwplot(matdev)
%title('stdev plot')
%caxis([15,126])

eval(['save img-',batch,' specmat filehist matdev']);
return

%%%%%%%%%%%below from spectrogram...clean up!%%%%%%%%
figure
newplot;

        args = {10*log10(abs(specmat)+eps)};
    % Axis labels
    xlbl = '';
    ylbl = '';

hndl = surf(args{:},'EdgeColor','none');

axis xy; axis tight;
colormap(jet);
caxis([3.5203,85.9716]);
% AZ = 0, EL = 90 is directly overhead and the default 2-D view.
view(0,90);

ylabel(ylbl);
xlabel(xlbl);
title(batch);
return;
