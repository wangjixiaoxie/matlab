%THIS SCRIPT IS ONLY FOR o57pu67 in 1/06

colordef black
figure

samp_fs=fs;
fs=44100;
rawsong=resample(rawsong,fs,samp_fs);
wn_start=3.2;
wn_end=3.35;



startsong=1;
endsong=6.5;


lowfiltval=75
highfiltval=20000

finlofiltval=150;
finhifiltval=12000;


%datinfiles={'6337_1.wav' '6337_2.wav'};
%aiffile='10055.aif';
%datoutfiles={'3667s1' '3667s2'};
%aif_outfiles={'10055.aif' '10055m10.aif' '10055p10.aif'}



datsong=rawsong(floor(fs*startsong):floor(fs*endsong));

%highpassfiltertoremovelow
[datsong]=bandpass(datsong,fs,lowfiltval,highfiltval,'butter');


normvec(1)=max(datsong);
normvec(2)=abs(min(datsong));
normfactor=max(normvec);

datsong=(1/(normfactor))*datsong;
pshift{1}=datsong;
%writeaiff_file
%aiffwrite(aiffile,datsong,fs,16);


%Write the aif files back to regular data%
%for(i=1:length(aif_outfiles))
 %   pshift{i}=aiffread(aif_outfiles{i});
%end

%for i=1:length(pshift)
    %pshift{i}=pshift{i}(1:262000)
    %[pshift{i}]=bandpass(pshift{i},fs,lowfiltval,highfiltval,'butter');
    %normvec(1)=max(pshift{i});
    %normvec(2)=abs(min(pshift{i}))

    %normfactor(i)=max(pshift{i});

    %ind{i}=find(pshift{i}==normfactor(i));
%end


%renormalize the data


figure;
numfiles=length(pshift)

%Examine the data
for(i=1:length(pshift))
    ax(i)=subplot(numfiles,1,i);
    [sm,sp,t,f]=evsmooth(pshift{i},fs,.05);
    imagesc(t,f,log(abs(sp)));syn;ylim([0,1e4]);
    smarray{i}=sm;
    spectarray{i}=sp;
    tarray{i}=t;
    farray{i}=f;
    %title();
    %axis([2 4 1 10000]);
    
end

linkaxes(ax,'x');

linkaxes(ax);

%Calculate the cross-correlation of the smoothed data
figure;
norm=smarray{1};
for(i=1:length(smarray))
       
        xcorrarray=xcorr(smarray{i},norm,5000);
        subplot(4,2,i)
        plot(-5000:5000,xcorrarray)
        offsetarray(i)=find(xcorrarray==max(xcorrarray))-5000;
end

%shift data back in time
for(i=1:length(pshift))
    %corpshift7980{i}=circshift(pshift7980{i}(:,1),-offsetarray{i});
    lng=length(pshift{i});

%instead cut out points
    tester=pshift{i};
    corpshift=tester(offsetarray(i)+1:lng);
%now add points at end%
    %if offset is 1, I want to just add the same point)
    sampvec=lng:-1:(lng-offsetarray(i)+1)
    pshiftbit=(tester(sampvec,1));
    corpshift=[corpshift ;pshiftbit];
    corpshifted{i}=corpshift;
    %check to make sure total power is constant
end

%check total power
figure
for(i=1:length(pshift))
    fftarray{i}=fft(corpshifted{i})
    fftfreqs=0:((44100/2)/(length(corpshifted{i})/2)):44100/2;
    subplot(4,2,i)
    absarray=abs(fftarray{i}).^2;
    plot(fftfreqs(1:100000), log(absarray(1:100000)))
end

%cut out initial points...6616 is a good blank spot
%for(i=1:length(pshift7980))
 %   corpshifted7980{i}=corpshifted7980{i}(6616:length(corpshifted7980{i}));
%end
%produce a white noise burst
start_spot=wn_start*fs;
end_spot=wn_end*fs;
randvec=2*rand(end_spot-start_spot+1,1)-1;
corpshifted{2}=corpshifted{1};
corpshifted{2}(start_spot:end_spot)=randvec;
corpshifted{3}=zeros(length(corpshifted{1}),1);
corpshifted{3}(start_spot:end_spot)=randvec;

numfiles=numfiles+2;
%produce a reverse bos stimulus
testvec=length(corpshifted{1}):-1:1;  
corpshifted{numfiles+1}=corpshifted{1}(testvec);

%THIS NORMALIZATION IS OUT OF ORDER%
%Then normalize
%for i=1:length(corpshifted)

%normvec(1)=max(corpshifted{i});
%normvec(2)=abs(min(corpshifted{i}));
%normfactor=max(normvec);


%corpshiftednorm{i}=corpshifted{i}.*(1/normfactor);
%end


%and weight the first 100 ms, 
%and the last 100 ms, with a

corpshiftednorm=corpshifted;

ln_endpoints=.1*fs
slope=1/ln_endpoints
inwtvector=0:slope:1
%fnwtvector=1:-slope:0;
begvector=1:1:length(inwtvector);
endvector=length(corpshiftednorm{1}):-1:length(corpshiftednorm{1})-length(inwtvector)+1;

for i=1:length(corpshiftednorm)
    corpshiftednormg{i}=corpshiftednorm{i}
    corpshiftednormg{i}(begvector)=corpshiftednorm{i}(begvector).*inwtvector'
    corpshiftednormg{i}(endvector)=corpshiftednorm{i}(endvector).*inwtvector'
end

%bandpass the stimuli again with a more stringent bandpass%
for i=1:length(corpshiftednormg)
    [corpshiftednormg{i}]=bandpass(corpshiftednormg{i},fs,finlofiltval,finhifiltval,'butter');
end

initbuff=zeros(2*fs,1);
finalbuff=zeros(fs,1);

for i=1:length(corpshiftednormg)
    [corpshiftednormg{i}]=[initbuff ;corpshiftednormg{i} ;finalbuff]
end

%Here do some stimulus checking, and create jiggled stimuli...


%FINALLY, WRITE TO WAV FILES*/
wavstring={'r78pu98.wav' 'r78pu98wn.wav' 'siwn.wav' 'r78pu98rev.wav'}
for i=1:length(corpshiftednormg)
    %corpshiftednormg{i}=[initbuff;corpshiftednormg{i};finalbuff];
    %resample(corpshiftednormg{i},44100,44150);
    wavwrite(corpshiftednormg{i}, 44100,16,wavstring{i});

end

