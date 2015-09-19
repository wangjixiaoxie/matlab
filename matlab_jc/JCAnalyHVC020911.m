
% .song.mat
    % Song data, sampled as 32kHz
% .spk
    % First column = cluster number (for each neuron, only a subset of
            % clusters are actual spikes - and the rest are just noise -
            % see e-mail from raghav)
    % Second column = spike time (seconds)
    
%     y50g89_081610160745.song.mat





% Undirected
load /cardinal/HVCAnalysedData/UDsaved021011.mat

% Directed
load /cardinal/HVCAnalysedData/Dsaved021011.mat

basedir='/cardinal/HVCAnalysedData/'
a=1;folders{a}=[basedir 'y50g89_081610_4/Directed'];clusters{a}=[1 3 4 5 6 7];prefixes{a}='Chan3_';
a=2;folders{a}=[basedir 'y50g89_081610_4/UnDirected'];clusters{a}=[1 3 4 5 6 7];prefixes{a}='Chan3_';
a=3;folders{a}=[basedir 'y50g89_081610_5/Directed'];clusters{a}=[1 2];prefixes{a}='Chan3_';
a=4;folders{a}=[basedir 'y50g89_081610_5/UnDirected'];clusters{a}=[1 2];prefixes{a}='Chan3_';
a=5;folders{a}=[basedir 'y50g89_081710_6/Directed'];clusters{a}=[3 4 5];prefixes{a}='Chan2_';
a=6;folders{a}=[basedir 'y50g89_081710_6/UnDirected'];clusters{a}=[3 4 5];prefixes{a}='Chan2_';
a=7;folders{a}=[basedir 'y50g89_081710_7/Directed'];clusters{a}=[2];prefixes{a}='Chan2_';
a=8;folders{a}=[basedir 'y50g89_081710_7/UnDirected'];clusters{a}=[2];prefixes{a}='Chan2_';
a=9;folders{a}=[basedir 'y50g89_081710_8/Directed'];clusters{a}=[1 6 9 12];prefixes{a}='Chan2_';
a=10;folders{a}=[basedir 'y50g89_081710_8/UnDirected'];clusters{a}=[1 6 9 12];prefixes{a}='Chan2_';
a=11;folders{a}=[basedir 'y50g89_081810_9/Directed'];clusters{a}=[3];prefixes{a}='Chan3_';
a=12;folders{a}=[basedir 'y50g89_081810_9/UnDirected'];clusters{a}=[3];prefixes{a}='Chan3_';
a=13;folders{a}=[basedir 'y50g89_082010_10/Directed'];clusters{a}=[4 7 9];prefixes{a}='Chan3_';
a=14;folders{a}=[basedir 'y50g89_082010_10/UnDirected'];clusters{a}=[4 7 9];prefixes{a}='Chan3_';
a=15;folders{a}=[basedir 'y50g89_082010_11/Directed'];clusters{a}=[1 6 9 10];prefixes{a}='Chan3_';
a=16;folders{a}=[basedir 'y50g89_082010_11/UnDirected'];clusters{a}=[1 6 9 10];prefixes{a}='Chan3_';
a=17;folders{a}=[basedir 'y50g89_082010_12/Directed'];clusters{a}=[2 4 7 9 12];prefixes{a}='Chan4_';
a=18;folders{a}=[basedir 'y50g89_082010_12/UnDirected'];clusters{a}=[2 4 7 9 12];prefixes{a}='Chan4_';
a=19;folders{a}=[basedir 'y50g89_082310_13/Directed'];clusters{a}=[2 5];prefixes{a}='Chan4_';
a=20;folders{a}=[basedir 'y50g89_082310_13/UnDirected'];clusters{a}=[2 5];prefixes{a}='Chan4_';
a=21;folders{a}=[basedir 'y50g89_082410_14/Directed'];clusters{a}=[5 7 10];prefixes{a}='Chan5_';
a=22;folders{a}=[basedir 'y50g89_082410_14/UnDirected'];clusters{a}=[5 7 10];prefixes{a}='Chan5_';
a=23;folders{a}=[basedir 'y50g89_101210_15/Directed'];clusters{a}=[1 3];prefixes{a}='Chan4_';
a=24;folders{a}=[basedir 'y50g89_101210_15/UnDirected'];clusters{a}=[1 3];prefixes{a}='Chan4_';



isDir=1; % for directed
[neuron]=HVCprocessing(folders,clusters,prefixes,isDir);
% only works for this bird!!!!
% neuron(folder)
%	 .song(song)
        % .Song --> entire song file
        % .noteons{syllable}(rendition)   % note onsets (minus 30ms) relative to song onset (30ms subtracted to allow pitch quant)
        % .noteoffs{syllable}(rendition)  % note offsets (minus 30ms) relative to song onset
        % .spiketimes{cluster}            % spike times (minus 30ms) relative to song onset (30ms subtracted to align with note onsets)
        % .notedata{syllable}.datt(rendition,:)   % raw data for note - used to calculate pitch contours
        % .pitchdata{syllable}.datt(rendition,:)  % pitch contours - first point is 16ms after note onset (512 points)


% PLOT RRRRRasters!
            pitchtimes=[1/8:1/8:1745/8]+16; % 16ms because of window size
            vals=1300;
            figure;hold on
            for thisfolder=1:11%length(neuron)
                goodclust=clusters{thisfolder*2};
                for thissyllable=1
                    for thissong=1:length(neuron(thisfolder).song)
                        plot(pitchtimes,neuron(thisfolder).song(thissong).pitchdata{thissyllable}.datt)
                        for thisclust=1:length(goodclust)
                            for thisrendition=1:length(neuron(thisfolder).song(thissong).noteons{thissyllable})
                            thesespikes=neuron(thisfolder).song(thissong).spiketimes{thisclust}...
                                -neuron(thisfolder).song(thissong).noteons{thissyllable}(thisrendition);
                            if ~isempty(thesespikes)
                                vals=vals-0.2;
                            plot(thesespikes,vals,'b.')
                            end
                            end
                        end
                    end
                end
            end



clearvars -except folders clusters prefixes
load /cardinal/HVCAnalysedData/UDsaved021011.mat
UDneuron=neuron;
clearvars -except UDneuron folders clusters prefixes
load /cardinal/HVCAnalysedData/Dsaved021011.mat
Dneuron=neuron;
clearvars -except UDneuron Dneuron folders clusters prefixes


wins=[-820:10:1420];
% Model 1: total number of spikes
    %%%%% This code does the original analyze - does total number of spikes in
    %%%%% a given window correlate with pitch in the harmonic stack?
    outvarsUD_model1=analyzeHVCallspikes(UDneuron,clusters,wins);
    outvarsD_model1=analyzeHVCallspikes(Dneuron,clusters,wins);
    %%%%%
    
% Model 2: spikes vs no spikes
    %%%%% Does the presence of any spikes in a given window correlate with
    %%%%% pitch in the harmonic stack?
    outvarsUD_model2=analyzeHVCbinaryspikes(UDneuron,clusters,wins);
    outvarsD_model2=analyzeHVCbinaryspikes(Dneuron,clusters,wins);
    %%%%%

% Model 3: total number of spikes (>0)
    %%%%% If we only consider renditions during which the neuron spiked
    %%%%% (i.e. ignoring renditions with no spikes), does the total number
    %%%%% of spikes in a given window correlate with pitch in the harmonic
    %%%%% stack?
    %%%%%
    outvarsUD_model3=analyzeHVCpositivespikes(UDneuron,clusters,wins);
    outvarsD_model3=analyzeHVCpositivespikes(Dneuron,clusters,wins);
    %%%%%



load /cardinal/HVCAnalysedData/forprint021111.mat

figure;hold on;
subplot(231);hold on;plot([-1000 1500],[0 0],'k')
plot([0 0],[-2 3],'k')
hold on;plot(wins,outvarsUD_model1.mnA-outvarsUD_model1.seA)
hold on;plot(wins,outvarsUD_model1.mnA+outvarsUD_model1.seA)
hold on;plot(wins,outvarsUD_model1.mnB+outvarsUD_model1.seB,'g')
hold on;plot(wins,outvarsUD_model1.mnB-outvarsUD_model1.seB,'g')
hold on;plot(wins,outvarsUD_model1.mnC+outvarsUD_model1.seC,'r')
hold on;plot(wins,outvarsUD_model1.mnC-outvarsUD_model1.seC,'r')
ylim([-2 3])
subplot(232);hold on;plot([-1000 1500],[0 0],'k')
plot([0 0],[-2 3],'k')
hold on;plot(wins,outvarsUD_model2.mnA-outvarsUD_model2.seA)
hold on;plot(wins,outvarsUD_model2.mnA+outvarsUD_model2.seA)
hold on;plot(wins,outvarsUD_model2.mnB+outvarsUD_model2.seB,'g')
hold on;plot(wins,outvarsUD_model2.mnB-outvarsUD_model2.seB,'g')
hold on;plot(wins,outvarsUD_model2.mnC+outvarsUD_model2.seC,'r')
hold on;plot(wins,outvarsUD_model2.mnC-outvarsUD_model2.seC,'r')
ylim([-2 3])
subplot(233);hold on;plot([-1000 1500],[0 0],'k')
plot([0 0],[-2 3],'k')
hold on;plot(wins,outvarsUD_model3.mnA)
hold on;plot(wins,outvarsUD_model3.mnB,'g')
hold on;plot(wins,outvarsUD_model3.mnC,'r')
ylim([-2 3])



subplot(234);hold on;plot([-1000 1500],[0 0],'k')
plot([0 0],[-2 3],'k')
hold on;plot(wins,outvarsD_model1.mnA-outvarsD_model1.seA)
hold on;plot(wins,outvarsD_model1.mnA+outvarsD_model1.seA)
hold on;plot(wins,outvarsD_model1.mnB+outvarsD_model1.seB,'g')
hold on;plot(wins,outvarsD_model1.mnB-outvarsD_model1.seB,'g')
hold on;plot(wins,outvarsD_model1.mnC+outvarsD_model1.seC,'r')
hold on;plot(wins,outvarsD_model1.mnC-outvarsD_model1.seC,'r')
ylim([-2 3])
subplot(235);hold on;plot([-1000 1500],[0 0],'k')
plot([0 0],[-2 3],'k')
hold on;plot(wins,outvarsD_model2.mnA-outvarsD_model2.seA)
hold on;plot(wins,outvarsD_model2.mnA+outvarsD_model2.seA)
hold on;plot(wins,outvarsD_model2.mnB+outvarsD_model2.seB,'g')
hold on;plot(wins,outvarsD_model2.mnB-outvarsD_model2.seB,'g')
hold on;plot(wins,outvarsD_model2.mnC+outvarsD_model2.seC,'r')
hold on;plot(wins,outvarsD_model2.mnC-outvarsD_model2.seC,'r')
ylim([-2 3])
subplot(236);hold on;plot([-1000 1500],[0 0],'k')
plot([0 0],[-2 3],'k')
hold on;plot(wins,outvarsD_model3.mnA)
hold on;plot(wins,outvarsD_model3.mnB,'g')
hold on;plot(wins,outvarsD_model3.mnC,'r')
ylim([-2 3])