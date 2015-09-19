function [specimg,matdev]=specspread1(batch)
%[specimg,devmat]=specspread1(batch)

%using lag=0. assuming Fs=44100 and w type file for evsoundin.
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

[s,f,t,p]=spectrogram(sounddata,512,300,[],44100);
clear s f t sounddata; pack;
maxfreq=fix(size(p,1)/2);
disp(['this file is ',num2str(size(p,2)),' fft segments long.'])
specmat(1:maxfreq+1,1:maxfreq+1,prevz+size(p,2))=zeros;
if filecount==1
	matdev(1+maxfreq,1+maxfreq)=0;
end

for z=1:size(p,2)
	%%progress report:
	if fix(z/500)==z/500	
		disp(['       ',num2str(z),'  ',num2str(prevz),'  ',num2str(size(specmat,2)),'  ',num2str(size(specmat,1))])
	end
	
	specmat(1,2:end,prevz+z)=p(1:maxfreq,z); %specmat(1,:,"z")is low-to-high freq spectr. unchanged.
	specmat(1,1,prevz+z)=1;
	specmat(2:end,1,prevz+z)=p(1:maxfreq,z); %specmat first col is scalars for freq bins = power in that band.

	for i=2:size(specmat,1)
		if z+lag > size(p,2) | z+lag <= 0
			specmat(2:end,2:end,prevz+z)=0;break; %danger!! zeros disrupt avg?
		end
	
	%	toadd=p(1:maxfreq,z+lag)'*specmat(i,1,prevz+z);
	%	if filecount~=1
	%		avgest=mean(specmat(i,2:maxfreq+1,1:prevz+z),3);
	%		matdev(i,2:end)=matdev(i,2:end)+(toadd-avgest).^2;
			
	%		else
	%		varz=z;
	%	end
		
	%	specmat(i,2:end,prevz+z)=toadd;
		specmat(i,2:end,prevz+z)=p(1:maxfreq,z+lag)'*specmat(i,1,prevz+z);
		%obsolete. wrong anyway: matdev(i,2:end)=matdev(i,2:end)+((p(1:maxfreq,z+lag)'*specmat(i,1,z))-(mean(specmat(i,2:end,1:prevz+z),3))).^2; %in specspread1, matdev was almost always about .88 to .98 times the real std
	end
end
filehist(filecount,:)=[filenm,' starts after ',prevz,', goes till ',prevz+z];
prevz=prevz+z;

end %%%end whileloop.. indent..%%%%
%specmat is built. average for output display

% bring back for normal output:
%specimg=mean(specmat,3); %now need to throw out the avg spec (row1) (and flip?)
%specimg=specimg(2:end,2:end);
specimg=specmat(2:end,2:end,:);

%specimg=flipud(specimg(2:end,:));

%DO STATS ON SPECMAT!! 
%clear specmat;

%eval(['save img-',batch,' specimg filehist']);

%non on-the-fly matdev. works great, is the same as using "std" function.
%matdev=specmat(:,:,1).*0;
%for z=1:size(specmat,3)
%	for i=2:size(specmat,1)
%	matdev(i,2:end)=matdev(i,2:end)+((p(1:maxfreq,z+lag)'*specmat(i,1,z))-(mean(specmat(i,2:end,:),3))).^2;
%	end
%end


%bring back to get a matdev again%%%%%%%%%%%%%%%%%%%%%
matdev=0;
%matdev=sqrt(matdev(2:end,2:end)./(prevz-varz-1));
%for CV instead of stdev add the next line:
%matdev=matdev./specmat; %this line to use matdev to report CV.
%%%%%%%%%%%%%%%%%%%

return

%%%%%%%%%%%below from spectrogram...clean up!%%%%%%%%
figure
newplot;

        args = {10*log10(abs(specimg)+eps)};
    % Axis labels
    xlbl = '';
    ylbl = '';

hndl = surf(args{:},'EdgeColor','none');

axis xy; axis tight;
colormap(jet);

% AZ = 0, EL = 90 is directly overhead and the default 2-D view.
view(0,90);

ylabel(ylbl);
xlabel(xlbl);
title(batch);
return;
