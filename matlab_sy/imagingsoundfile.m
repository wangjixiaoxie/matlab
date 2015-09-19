%THIS SCRIPT IS ONLY FOR o57pu67 in 1/06

colordef black

figure

rawsong{1}=g88s;
rawsong{2}=pu81pk52s;

for i=1:2
    rawsong{i}=rawsong{i}-mean(rawsong{i});
end


start_times=[3.45 2.09];
song_lng=7;
wavstring={'g88s' 'pu81pk52s'};

samp_fs=44150;
fs=44100;
rawsong{1}=resample(rawsong{1},fs,samp_fs);
rawsong{2}=resample(rawsong{2},fs,samp_fs);

for i=1:2
    datsong{i}=rawsong{i}(floor(fs*start_times(1)):floor(fs*(start_times(1)+8)));
end
%highpassfiltertoremovelow

initbuff=zeros(2*fs,1);
finalbuff=zeros(fs,1);



ln_endpoints=.1*fs
slope=1/ln_endpoints
inwtvector=0:slope:1
%fnwtvector=1:-slope:0;
begvector=1:1:length(inwtvector);
endvector=length(datsong{1}):-1:length(datsong{1})-length(inwtvector)+1;

for i=1:2
    datsong{i}=datsong{i}
    datsong{i}(begvector)=datsong{i}(begvector).*inwtvector'
    datsong{i}(endvector)=datsong{i}(endvector).*inwtvector'
end

for i=1:length(datsong)
    [datsong{i}]=[initbuff ;datsong{i} ;finalbuff]
    datsong{i}=datsong{i}/(1.01*max(abs(datsong{i})));
    wavwrite(datsong{i},44100,16,wavstring{i});
end


