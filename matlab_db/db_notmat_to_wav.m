function [] = db_notmat_to_wav(batchfile, foldername)

mkdir(['../' foldername])

fid = fopen(batchfile, 'r');
tline = fgetl(fid);

while ischar(tline)
	db_normalize_song_waveform(tline,1);
	movefile([tline(1:end-4) 'wav'], ['../' foldername])
	tline = fgetl(fid);
end
