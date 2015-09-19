function mov=specmovie(anchor,batchfile,pre,post,tag)
% mov=specmovie(anchor,batchfile,pre,post,tag)
% makes movie where the frames are spectrograms of songs.
% mov is the movie object's handle
% anchor is a label found in the 'labels' variable in the .not.mat files for the songs in batchfile
% pre and post are # of seconds before and after the anchor to be included in the movie frames
% batchfile can be wavs or .not.mat files corresponding to wavs
% tag is a string that identifies the movie-- printed in each frame w./ frame#.

index=0; %the row-index of the output matrix
offst=0;
filenumb=0;
frame=0;
movieh=avifile([batchfile,'-',num2str(pre),'-',num2str(post),'movie.avi']);
fid=fopen(batchfile);
while 1
soundfile=fgetl(fid);
if ~ischar(soundfile), break, end
if frame~=0
	fprintf(1,'\r')
end
if soundfile(length(soundfile)-7:length(soundfile))=='.not.mat'
	if ~exist([soundfile])
		disp([soundfile,' not found.'])
		fclose(fid);
		return
	end
	load([soundfile]);
	disp(['computing ' soundfile])
	soundfile=soundfile(1:length(soundfile)-8);
	else
	if ~exist([soundfile,'.not.mat'])
		disp([soundfile,'.not.mat not found.'])
		fclose(fid);
		return
	end
	load([soundfile,'.not.mat']);disp(['computing ' soundfile '.not.mat'])
end
[song, Fs]=evsoundin('',soundfile,'w');
filenumb=filenumb+1;
anchors=find(labels==anchor);

if numel(anchors)==0
continue
end

offst=0;		%is that right? should be 1?
if floor((Fs/1000)*(onsets(anchors(1))-(1000*pre)))<=0
	offst=1+abs(floor((Fs/1000)*(onsets(anchors(1))-(1000*pre))));
	fprintf(1, '(file is short)')			% WARNING THERE MAY BE DIFFS B/W THESE VALUES AND WHAT IS USED IN 
end							% MAKING THE SPECGRAMS

song=[zeros(offst,1); song];
onsets=onsets+offst/(Fs/1000);
offsets=offsets+offst/(Fs/1000);


for i=1:length(anchors)
%if i~=1, fclose(fid),return,end
	frame=frame+1;
	fprintf(1,['.. ',num2str(frame)])
	if floor((Fs/1000)*(onsets(anchors(i))-(1000*pre)))<0
		fprintf(1, '(WHAAAAAT?? still short)')
		offst=abs(floor((Fs/1000)*(onsets(anchors(i))-(1000*pre))));	%DELETE ALL B/W IF AND ELSE
		song=[zeros(offst,1); song];					%IF 'WHAAAAT??' NEVER APPEARS.
		strt=1;
		stp=offst+ceil((Fs/1000)*(offsets(anchors(i))+(1000*post)));
	else
		strt=ceil((Fs/1000)*(onsets(anchors(i))-(1000*pre)));
		stp=ceil((Fs/1000)*(offsets(anchors(i))+(1000*post)));
	end
	if ceil((Fs/1000)*(offsets(anchors(i))+(1000*post)))>length(song)	% maybe this doesn't need to be in an if loop--1st zeros arg may be 0 almost always
		song=[song; zeros(ceil((Fs/1000)*(offsets(anchors(i))+(1000*post)))-length(song),1)];
		fprintf(1, '(long)')
	end
	%Find axis labels: start, anchor, end, anchor +1, etc........
	%fprintf(1,[':',num2str(strt),' ',num2str(stp)])
	framedat=song(strt:stp);
	fig=figure('visible','off');
	specgram(framedat,512,Fs,[],384);
	caxis(((2^8)-1)*[0, .4]);
	text(0,10000,[tag,'- ',num2str(frame)],'fontsize',16);
	movieh=addframe(movieh,fig);
%	eval(['print -djpeg movieframe',num2str(frame)]);
	close(fig);
end
end %for while loop - (batchfile fgetl)
fclose(fid);
movieh=close(movieh);
figure;
mov=aviread([batchfile,'-',num2str(pre),'-',num2str(post),'movie.avi']);
movie(mov)
