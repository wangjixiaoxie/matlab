function make_mahal_stim(mahal_strct,name)
%
%function make_mahal_stim(mahal_strct,name)
%outputs a krankfile with the (32khz, binary sound file)
%
%mahal_strct is a structure that is in acitve in your workspace (for now)
%
%name is the name of the stimulus file to be created
%

num_syls = size(mahal_strct);
num_syls = num_syls(2);
fs = mahal_strct{1}.FS;
z_pad = zeros(1,round(60*(fs/1000)));
z_pads = 60;
sng = z_pad;
labels = [];
cnt = 0;
for i = 1:num_syls
    for j = 1:10
        cnt = cnt+1;
        crntsyl = mahal_strct{i}.syl{j}.snd;
        crntlng = length(crntsyl)*((fs/1000)^-1);
        sng = [sng z_pad crntsyl'];
        labels(i) = mahal_strct{i}.label;
        if i ==1
            onsets(cnt) = z_pads;
            offsets(cnt) = onsets(cnt)+crntlng;
        else
            onsets(cnt) = offsets(cnt-1)+z_pads;
            offsets(cnt) = onsets(cnt)+crntlng;
        end
    end
end
labels  = char(labels);
sng = [sng z_pad];
sng = resample(sng,32e3,fs);
namesng = [char(name) '.mahal_syls'];
%krankwrite(char(namesng),sng);
make_obs(char(namesng),sng);
namenotmat = [char(namesng) '.not.mat'];
save(char(namenotmat), 'labels','onsets','offsets');    