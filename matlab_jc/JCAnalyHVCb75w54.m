
% .song.mat
    % Song data, sampled as 32kHz
% .spk
    % First column = cluster number (for each neuron, only a subset of
            % clusters are actual spikes - and the rest are just noise -
            % see e-mail from raghav)
    % Second column = spike time (seconds)
    
%     y50g89_081610160745.song.mat





% Undirected


% Directed


basedir='/cardinal/HVCAnalysedData/'
a=1;folders{a}=[basedir 'b75w54_032010_1/Directed'];clusters{a}=[2 4 8 10];prefixes{a}='Chan4_';
a=2;folders{a}=[basedir 'b75w54_032010_1/UnDirected'];clusters{a}=[2 4 8 10];prefixes{a}='Chan4_';




isDirected=1; % for directed
[Dneuron]=HVCprocessing_b75w54(folders,clusters,prefixes,isDirected);
isDirected=0; % for directed
[UDneuron]=HVCprocessing_b75w54(folders,clusters,prefixes,isDirected);

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
            vals=1250;
            figure;hold on
            for thisfolder=[1 3]%length(neuron)
                goodclust=clusters{thisfolder*2};
                for thissyllable=1
                    for thissong=1:length(neuron(thisfolder).song)
                        plot(pitchtimes,neuron(thisfolder).song(thissong).pitchdata{thissyllable}.datt)
                        for thisclust=goodclust
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
load /cardinal/HVCAnalysedData/UDsaved_b75w54.mat
clearvars -except UDneuron folders clusters prefixes
load /cardinal/HVCAnalysedData/Dsaved_b75w54.mat
clearvars -except UDneuron Dneuron folders clusters prefixes

for i=1:length(UDneuron)
    ALLneuron(i*2)=UDneuron(i);
    ALLneuron(i*2-1)=Dneuron(i);
end
harmwin=[50 70;130 160;100 160];

wins=[-820:10:1420];
% Model 1: total number of spikes
    %%%%% This code does the original analyze - does total number of spikes in
    %%%%% a given window correlate with pitch in the harmonic stack?
    [outvarsUD_model1,outdataUD_model1]=analyzeHVCallspikes2(ALLneuron,clusters,wins,harmwin,[2],1);
    [outvarsD_model1,outdataD_model1]=analyzeHVCallspikes2(ALLneuron,clusters,wins,harmwin,[1],1);    
    