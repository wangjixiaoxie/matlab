bos='D:\tim\w90o24\stim\9024.wav'
rev=       'D:\tim\w90o24\stim\9024rev.wav'
numbos=3;
ff=dir('*.cbin');
song_chan=1;
tet_chans=2:5;
i=9;
while(numbos==3)
    fn=ff(end-i).name
    rd=readrecf(fn)
    if (strcmp(rd.outfile,rev))
        break;  
    
    else
       fn=ff(end-i).name
        i=i+1;
    end
    
end
[data,fs,spkindnew,spkampnew]=tetanaldc(fn,-1500,song_chan,tet_chans);
spkind=spkindnew;
spkamp=spkampnew;
pltdat;

pltscat;
