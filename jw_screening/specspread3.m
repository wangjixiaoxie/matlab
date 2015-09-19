function [spectra]=specspread3(batch,savename)
%[spectra]=specspread3(batch,savename)
% specspread3 is a compromise between 1(saves full distributions) and 2(gets means and std only).
% in specspread3, nothing is computed-- but all spectral segments are saved into a matrix so that
% the specmat for any lag can be computed from it. It basically saves spectrograms for each file in
% the batchfile. .mat has batch name AND savename in the filename.
%NB: why is specspread2 slower than ss1 or 3? Variable mem allocation.. .fixed in ss1 I guess.
%oldnotes:
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
	[sounddata, Fs] = evsoundin('',filenm,'w');
	sounddata = highpass(sounddata,200,Fs);
else
	disp(['file ', filenm,' doesn''t exist on this path?']); continue;
end

[s,f,t,p]=spectrogram(sounddata,512,300,[],44100); %sounddata,512,300,[],44100 means 212samples per segment.
maxfreq=fix(size(p,1)/2);
disp(['this file is ',num2str(size(p,2)),' fft segments long.'])
spectra(1,1:maxfreq,prevz+size(p,2))=zeros;

for z=1:size(p,2)
	%%progress report:
	if fix(z/500)==z/500	
		disp(['       ',num2str(z),'  ',num2str(prevz)])
	end
	
	if z+prevz==1;
	%specmat(maxfreq+1,maxfreq+1)=0;
	%matdev(maxfreq+1,maxfreq+1)=0;
	end
	spectra(1,1:maxfreq,z+prevz)=p(1:maxfreq,z)';
	%specmat(1,2:end)=specmat(1,2:end)+p(1:maxfreq,z)'; 
	%specmat(2:end,1)=p(1:maxfreq,z);
	
end
filehist(filecount,:)={filenm,'starts after',num2str(prevz),'goes till',num2str(prevz+z)};
prevz=prevz+z;

end %%%end whileloop.. indent..%%%%

%specmat=(specmat./prevz);
%avgspec=specmat(1,2:end);
%specmat=specmat(2:end,2:end);
%specimg=flipud(specimg(2:end,:));
%matdev=sqrt(matdev(2:end,2:end)./(prevz-varz-1)); %to use matdev to report stdev, comment out next line. Checked with highflat.. this is a good std estimate.
%matdev=matdev./specmat; %this line to use matdev to report CV.
%%specmat(5000:100:7000)
%%matdev(5000:100:7000)
%%jwplot(matdev)
%%title('stdev plot')
%%caxis([15,126])

%eval(['save sssegments-',batch,'-',savename,' spectra filehist']);
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
