%eval('!ls -1 *.not.mat > notmat_list.txt');
MAX_NOTES = 10;
srch_str = 'b';
av_dur = [-150,250]; % duration of average in ms

fid = fopen('batchnotmat','r');
fn  = fgetl(fid);
fclose(fid);

Fs = 44100.0;Flw = 300.0;Fhi = 8000.0; %in Hz
nfft = 512; spectwin = nfft; noverlap = floor(nfft*0.8);

basename = fn(1:end-8);
[rsong,Fs]=evsoundin('',basename,'ebin0r');
load([basename,'.not.mat']);
filtsong = bandpass(rsong,Fs,Flw,Fhi,'hanningfir');
[sp,f,t] = specgram(filtsong,nfft,Fs,spectwin,noverlap);
dt = t(2)-t(1);
dat_dur = ceil(av_dur*1e-3/dt);

notes = zeros([floor(nfft/2)+1,diff(dat_dur)+1,MAX_NOTES]);
Nnotes = zeros([MAX_NOTES,1]);

fid = fopen('batchnotmat','r');
while (1)
	fn = fgetl(fid);
	if (~ischar(fn))
		break;
	end
	basename = fn(1:end-8);
	disp(basename);
	if (exist(basename,'file')&exist(fn,'file'))
		% load up labels and spectrogram
		[rsong,Fs]=evsoundin('',basename,'obs0r');
		load([basename,'.not.mat']);

		%[sp,nfft,spwin,nov,tmn,tmx,fmn,fmx]=read_spect(fn2);
		filtsong = bandpass(rsong,Fs,Flw,Fhi,'hanningfir');
		[sp,f,t] = specgram(filtsong,nfft,Fs,spectwin,noverlap);
		cmd = ['save ',basename,'.espect f t nfft Fs spectwin noverlap'];eval(cmd);

		% change onsets to indicies
		dt = t(2)-t(1);
		onsets  = ceil(onsets*1e-3/dt);

		%p = find(sp==max(max(sp)));sp(p) = min(min(sp));
		%p = find(sp <= 5);sp2=sp;sp2(p) = 2;

		s_pos = findstr(labels,srch_str);
		for ii = 1:length(s_pos)
			ind = s_pos(ii);nind = 1;
			sp_inds = onsets(ind) + [dat_dur(1):dat_dur(2)];
			notes(:,:,nind) = notes(:,:,nind)+abs(sp(:,sp_inds));
			Nnotes(nind) = Nnotes(nind) + 1;
			ind = ind + 1;nind = nind + 1;
		end
	end
end
fclose(fid);
