%%%% Only use /swift 5, /swift6, /swift7, /swift8, /swift9

%%% For general purposes
        % Make a directory of the song files
        dirf('*.cbin','batch');

        % To select only the catch trials (without WN)
        findcatch('batch') % this creates 'batch.catch'

        % To select only a certain proportion (e.g. 0.1) of the files
        randsamp('batch',0.1)

        % Open the GUI, segment, and label song - this creates *.not.mat files
        evsonganaly
            % Default segmentation is normally okay.
        % If many of the files are noisy, you may want to throw away the bad files
     cleandir4('batch',100000,500,5,5);  
            %   cleandir4('batch',10000,500,6,10);  
                %%% creates batch.keep and batch.dcrd
                %%% use evsonganaly to make sure batch.dcrd is garbage
            %   mk_rmdata('batch.dcrd',1)
            %   !csh yes | ./rmdata
                %%% do it for .tmp,.rec as well - by going into batch.dcrd and
                %%% replace all .cbin with .tmp

        % Make a directory of the *.not.mat files
        % RENAME
        % RENAME
        % RENAME
        % RENAME
        % RENAME
        % RENAME files before executing to avoid overwriting
        
     dirf('*.cbin.not.mat','batchnotes')
     fvalsOFF0415=findwnoteJC('batchnotes','a','','',0,[2000 2700],8500,1,'obs0',1);
     pitchOFF0415=pitchcontour(fvalsOFF0415,2000,3000);
     figure;plot(pitchOFF0415)
     wind=300:1:600;
     mnpitchOFF0415=mean(pitchOFF0415(wind,:));
     figure;plot(mnpitchOFF0415)
     figure;hist(mnpitchOFF0415)
     prctile(mnpitchOFF0415,70)
     
     %To look at the mean pitch
    
     figure;hold on;plot(mean(NAME'),'b');plot(mean(NAME'),'r')
     
     vals0407=evtaf_freq('batch.keep.rand',[2000,3000],'a',128,'obs0',0,0);
     FFvals0407=vals0407(:,2);
     prctile(FFvals0407,70)
        % Look at the average spectral structure of a labeled note - choose FF
            % range around a powerful harmonic
        edit tafsimTD
        
        edit mkwn
   % autolabel
   load Templates040611.mat
    mk_tempf('batch.catch',templaA,2,'obs0');    
    get_trigt2('batch.catch',cntrngA,0.2,128,1,1);   
    label_trigs('batch.catch','a','obs0',10000,1,1,5,30);
    
    % analyze
     fvals0407wnon=findwnoteJC('batchnotes','a','','',0,[2000 2700],8500,1,'obs0',1);
     pitch0407wnon=pitchcontour(fvals0407wnon,2000,3000);
  
     
%make wn file
%length of stimulus.
t=.05
fs=44100
npts=fs*t
rawdataflt=rand(npts,1)*2-1


    %normalize
    mxvl=max(abs(rawdataflt))
    rawdatafltnrm=rawdataflt/(mxvl*1.0001);  
    wavwrite(rawdatafltnrm,44100,16,'filename.wav')

    % times=timing3(fvalsOFF0415);
    % figure;plot(times,mnpitchOFF0415)
     
    