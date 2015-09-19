% Only use /swift 5, /swift6, /swift7, /swift8, /swift9

%%% For lesion analysis
LesionAnalysisRotation.m



%%% For general purposes
        % Make a directory of the song files
        dirf('*.cbin','batch') % dirf('*.wav','batch')

        % To select only the catch trials (without WN)
        findcatch('batch') % this creates 'batch.catch'

        % To select only a certain proportion (e.g. 0.1) of the files
        randsamp('batch',0.5)

        % Open the GUI, segment, and label song - this creates *.not.mat files
        evsonganaly
            % Default segmentation is normally okay.
        % If many of the files are noisy, you may want to throw away the bad files
            %   cleandir4('batch',10000,500,6,10);  
                %% creates batch.keep and batch.dcrd
                %% use evsonganaly to make sure batch.dcrd is garbage
            %   mk_rmdata('batch.dcrd',1)
            %   !csh yes | ./rmdata
                %%% do it for .tmp,.rec as well - by going into batch.dcrd and
                %%% replace all .cbin with .tmp

        % Make a directory of the *.not.mat files
        dirf('*.cbin.not.mat','batchnotes')


        % Analyze data - takes about 40 seconds per 100 syllables
        [fvalsPOST,tvalsPOST,pitchcurvesPOST]=summary_stats2010('batchnotes',1,[2800 3800],'a');
        % [2000 3000] is good for long stack notes 
        % [2800 3800] is good for short stack notes
        figure;plot(pitchcurvesPOST)
        figure;plot(mean(pitchcurvesPOST'))
        figure;hold on;subplot(122);hold on;plot(tvalsPOST,mean(pitchcurvesPOST(1000:1300,:)),'*','Color','r')
        subplot(121);plot(tvals,mean(pitchcurves(1000:1300,:)),'*')



        % Look at the average spectral structure of a labeled note - choose FF
            % range around a powerful harmonic
        edit tafsimOA