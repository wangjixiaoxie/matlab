for i=1:10
    ind=find(spkind(:,3)==i)
    if(isempty(ind)==0)
        ifn(i)=spkind(ind(1),2);
        fn=fnm{ifn(i)}
        [data,fs,spkindnew,spkampnew]=tetanaldc(fn,-2000,song_chan,tet_chans);
        datalengthcheck(i)=length(data(:,1));
        rd=readrecf(fn);
        stimlistcheck{i}=rd.outfile;

    end
end