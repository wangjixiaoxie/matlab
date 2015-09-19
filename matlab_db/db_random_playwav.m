% This script will play the .wav files in your current directory
% All the .wav files need to be listed in a .txt file called 
% 'song_files.txt'

% This opens the file 'song_files.txt'
fid = fopen('song_files.txt','r');

% Gets the first line of 'song_files.txt'
tline = fgetl(fid);
i = 1;
songlist = {};

% Loop goes through each line in song_files and saves
% it to the variable songlist
while ischar(tline)
	songlist{i} = tline;
	tline = fgetl(fid);
	i = i+1;
end

% This is the number of playback songs
number_of_songs = 30;

% This selects a random subset of the songs from the songlist
playlist = randi(length(songlist), [number_of_songs, 1]);

% Plays each song on playlist and then pauses for 20 seconds
% between songs. Use CRTL-C to escape.
for i = 1:length(playlist)
	current_song = wavread(songlist{playlist(i)});
	sound(current_song, 32000);
	display(['Playing ' songlist{playlist(i)}])
	display(['You have ' num2str(length(playlist)-i) ' songs to go.'])
	pause(10);
	display('Pause between songs');
	display(' ');
	pause(10);
end
	

