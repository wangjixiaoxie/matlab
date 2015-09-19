function [songnumber,songposition,notetime]=jctiming
% songnumber - this note was sung in the Nth song in the batch file
% songposition - this note was the Nth note of that type sung in the song
% notetime - this note was sung at N seconds after midnight

% First load a batch of note files, then load a batch of rec files

% Calculate times of song onsets
[labels,onsets,offsets]=notestats; % takes note files
varargout=evloadfile; % takes rec files
for i=1:length(varargout)
    songonset(i)=songtime(varargout(i).fname);
end

% Calculate the number of the song that each note appears in
indstart=find(labels=='/');
indnotes=find(labels=='a');
for i=1:length(indnotes)
    songnumber(i)=max(find(indstart<indnotes(i)));
end

% Calculate position of each note in its song
songposition(1)=1;
for i=2:length(songnumber)
    if songnumber(i)==songnumber(i-1)
        songposition(i)=songposition(i-1)+1;
    else
        songposition(i)=1;
    end
end

% Calculate the time of day that the note was sung
for i=1:length(indnotes)
    notetime(i)=songonset(songnumber(i))+onsets(indnotes(i))/1000;
end


g=7;