% LT modified 8/5 to incorporate into the "lt_compare_hits_to_misses"
% script.  specifically changed save structure since have 'b'(hits) and
% 'N' (misses).   removed the ende of this original script looing at
% individual slices.  

% LT 7/27/13 - modified from lt_db_make_template, which takes an average
% spec over all instances of a syllab eto make a template.  Here I will
% also make a cell of individual spectrograms, then for each count the
% number of peaks using findpeaks.



%% Choosing batch file and syllables; making variables

clear make_temp

% make_temp.parameters.batch_name = input('What is the name of the batch file?  ', 's');
% make_temp.parameters.birdname = input('What is the name of the bird? (bk90bk46)  ', 's');
batch_yes_or_no=input('use batch.catch.keep? (y or n): ', 's');
if batch_yes_or_no=='y';
    make_temp.parameters.batch_name = 'batch.catch.keep';
else
    make_temp.parameters.batch_name=input('what is the batch name? (no quote marks) ','s');
end    
make_temp.parameters.birdname = 'pu13bk43';
% % make_temp.parameters.syllables = input('What are the syllables? (no spaces)  ', 's');
make_temp.parameters.syllables='bN'

make_temp.spec = cell(1,length(make_temp.parameters.syllables));
make_temp.time_range = cell(1,length(make_temp.parameters.syllables));
make_temp.freq_range = cell(1, length(make_temp.parameters.syllables));
make_temp.template = cell(1, length(make_temp.parameters.syllables));
make_temp.parameters.slices = cell(1,length(make_temp.parameters.syllables));



% dir for saving
% suffix=fix(clock);
% mkdir(['lt_count_harmonics_results' num2str(suffix(1)) num2str(suffix(2)) num2str(suffix(3)) '_' num2str(suffix(4)) num2str(suffix(5))]);
% save_dir=[pwd '/lt_count_harmonics_results' num2str(suffix(1)) num2str(suffix(2)) num2str(suffix(3)) '_' num2str(suffix(4)) num2str(suffix(5))];


%% Added by LT, in order to look specifically at syallables followed or preceded by other specific syllables.  (e.g. only look at b if preceded by a)

% make_temp.parameters.prent=input('Enter syllables that must come before (blank if none): ', 's');
% make_temp.parameters.postnt=input('Enter syllables that must come after (blank if none): ', 's');

make_temp.parameters.prent='';
make_temp.parameters.postnt='';


%% Get AVN

% Makes two spectrograms around your syllable of interest. Used to choose
% the points for making your slice (template)
for i = 1:length(make_temp.parameters.syllables);
    [make_temp.spec{i},make_temp.time_range{i},make_temp.freq_range{i},~,make_temp.num_of_renditions(i)] = lt_get_avn(make_temp.parameters.batch_name,...
        make_temp.parameters.syllables(i),...
        0.2,0.2,make_temp.parameters.prent,make_temp.parameters.postnt,'obs0');
    
       
    figure(1); subplot(2,5,5*i-4); hold on
    imagesc(make_temp.time_range{i}*1000,make_temp.freq_range{i},log(make_temp.spec{i}));syn;ylim([0,1e4]);xlim([-150 150])
    title(['Syllable ' make_temp.parameters.syllables(i)]) % ', preceded by: ' make_temp.parameters.prent '; followed by: ' make_temp.parameters.postnt ])
    ylabel('Frequency (Hz)')
    xlabel('Time (ms)')
    
    figure
    imagesc(log(make_temp.spec{i}));syn;
    title(['Syllable ' make_temp.parameters.syllables(i) ', preceded by: ' make_temp.parameters.prent '; followed by: ' make_temp.parameters.postnt ])
    ylabel('Frequency (Hz)')
    xlabel('Points')    
end

disp('number of renditions [hit miss]')
disp(make_temp.num_of_renditions)

%% Get all individual spectrograms and put into a cell array.
% exactly ike get avn above except modified to keep all individual specs,
% instead of averaging.

clear make_temp.all_spectrograms

for i = 1:length(make_temp.parameters.syllables)
    [make_temp.all_spectrograms{i},make_temp.time_range{i},make_temp.freq_range{i}] = lt_get_individual_avn(make_temp.parameters.batch_name,...
        make_temp.parameters.syllables(i),...
        0.2,0.2,make_temp.parameters.prent,make_temp.parameters.postnt,'obs0');
end


%% Choose slices from the average spec.
% will use the slice coordinates not on the average, but on each individual spec, and then put all of those slices into a cell array.

clear peak_statistics
clear peaks


for i = 1:length(make_temp.parameters.syllables)
    % asks for slice points to construct template
    make_temp.parameters.slices{i} = input(['What are the slice coordinates for ' ...
        make_temp.parameters.syllables(i) '?\n(i.e. 60 or 70:82)  ']);
    for ind=1:size(make_temp.all_spectrograms{i},3)
        temporary_slice=mean(make_temp.all_spectrograms{i}(1:2:256,make_temp.parameters.slices{i},ind),2);
        if ind==1;
            size(temporary_slice) % should be 128 point vector
        end
        temporary_slice(1:6)=0;
        temporary_slice=temporary_slice./max(temporary_slice);
        make_temp.all_spectrogram_slices{i}(:,ind)=temporary_slice;
        clear temporary_slice
        
        
    end
% plot average slice

    make_temp.template{i} = mean(make_temp.spec{i}(1:2:256,make_temp.parameters.slices{i}),2);
    size(make_temp.template{i}) %should be 128 point vector
    make_temp.template{i}(1:6) = 0;
    make_temp.template{i} = make_temp.template{i}./max(make_temp.template{i});
    
    figure(1); subplot(2,5,5*i-3); hold on
    plot(make_temp.template{i}), title(['syllable ' make_temp.parameters.syllables(i)]); xlim([0 128])

    
    %plot all arbitrary temporary slice
        figure(1); subplot(2,5,5*i-2); hold on;  plot(make_temp.all_spectrogram_slices{i}), title(['syllable ' make_temp.parameters.syllables(i) ' a single instance.'])
%         saveas(figure(gcf), [save_dir '/all_spectrogram_slices_' make_temp.parameters.syllables(i)],'fig'); 
        xlim([0 120])
        
%         %plot 1 out of every 10 slices
%         figure, plot(make_temp.all_spectrogram_slices{i}(:,1:10:size(make_temp.all_spectrogram_slices{i},2))), title(['syllable ' make_temp.parameters.syllables(i) ' a single instance.'])
%         saveas(figure(gcf), [pwd '/all_spectrogram_slices_' make_temp.parameters.syllables(i) '_skipping10'],'fig')
% 
%         % plot 1 out of every 25 slices
%         figure, plot(make_temp.all_spectrogram_slices{i}(:,1:25:size(make_temp.all_spectrogram_slices{i},2))), title(['syllable ' make_temp.parameters.syllables(i) ' a single instance.'])
%         saveas(figure(gcf), [pwd '/all_spectrogram_slices_' make_temp.parameters.syllables(i) '_skipping25'],'fig')

        
        % find peaks in slices   
        peak_statistics.peak_height=0.12 % input('what is minimum peak height (range 0-1)? ');   
    for ind=1:size(make_temp.all_spectrograms{i},3);
    
    [peaks,locations]=findpeaks(make_temp.all_spectrogram_slices{i}(:,ind),'MINPEAKHEIGHT',peak_statistics.peak_height);
        make_temp.all_spectrogram_peakheights{i}{ind}=peaks;
        make_temp.all_spectrogram_peaklocations{i}{ind}=locations;
        make_temp.all_spectrogram_peakamounts{i}(ind)=length(peaks);
        
      
    end
            % getting rel amplitudes of 1st 3 harmonics
    %define harmonics (based on inspection of all slices overlayed in an open figure)
            harm= 'd' %input('what x-coordinates bound your harmonics? (enter "d" for default: [22 32; 47 59; 76 84])');
           if harm=='d';
               harm=[22 32; 47 59; 76 84];
           end
    for kk=1:3; %kk = harmonic #    
        for ind=1:size(make_temp.all_spectrograms{i},3);

            make_temp.all_spectrogram_ampl_first_three_harmonics{i}(kk,ind)=max(make_temp.all_spectrogram_slices{i}(harm(kk,1):harm(kk,2),ind));
        end
        
            %mean and std amplitude for each harmonic
            make_temp.all_spectrogram_mean_ampl_three_harmonics{i}(kk)=mean(make_temp.all_spectrogram_ampl_first_three_harmonics{i}(kk,:));
            make_temp.all_spectrogram_STD_ampl_three_harmonics{i}(kk)=std(make_temp.all_spectrogram_ampl_first_three_harmonics{i}(kk,:));

    end
    
        % for each syllable, calculate spectral mean using the three harmonics (i.e. freq that is
        % midpoint of the three harmonics, weighted by their amplitudes.)
        % (e.g. gives range from 1-3 depending on which harmonic dominates.
        make_temp.all_spectrogram_spectral_means{i}=nan(1,size(make_temp.all_spectrograms{i},3));
        for ind=1:size(make_temp.all_spectrograms{i},3);
             make_temp.all_spectrogram_spectral_means{i}(ind)=(sum((make_temp.all_spectrogram_ampl_first_three_harmonics{i}(:,ind).*[1;2;3]),1))./(sum(make_temp.all_spectrogram_ampl_first_three_harmonics{i}(:,ind),1));
        end
        
        % get mean and std over trials.  
        make_temp.mean_over_renditions_spectral_mean{i}=mean(make_temp.all_spectrogram_spectral_means{i}(:));
        make_temp.STD_over_renditions_spectral_mean{i}=std(make_temp.all_spectrogram_spectral_means{i}(:))

        % plot mean and std for each harmonic (across syllables), along
        % with spectram mean as a vertical line.
        figure(1); subplot(2,5,5*i-1); hold on
        errorbar([1 2 3], make_temp.all_spectrogram_mean_ampl_three_harmonics{i},make_temp.all_spectrogram_STD_ampl_three_harmonics{i});xlim([0.9 3.1])
        line1=line([make_temp.mean_over_renditions_spectral_mean{i} make_temp.mean_over_renditions_spectral_mean{i}],[0.01 1.2]);
        line2=line([make_temp.mean_over_renditions_spectral_mean{i}+make_temp.STD_over_renditions_spectral_mean{i} make_temp.mean_over_renditions_spectral_mean{i}+make_temp.STD_over_renditions_spectral_mean{i}],[0.01 1.2]), set(line2,'LineStyle','--');
        line3=line([make_temp.mean_over_renditions_spectral_mean{i}-make_temp.STD_over_renditions_spectral_mean{i} make_temp.mean_over_renditions_spectral_mean{i}-make_temp.STD_over_renditions_spectral_mean{i}],[0.01 1.2]), set(line3,'LineStyle','--');

        title({['mean amplitude of harmonics. syllable ' make_temp.parameters.syllables(i)],['vertical lines indicate spectral mean +/ STD']})
%         saveas(figure(gcf),[save_dir '/harmonics_mean_and_std_' make_temp.parameters.syllables(i)], 'fig') 
       
    
        %plot histogram of peak amounts
        figure(1); subplot(2,5,5*i); hold on
        hist(make_temp.all_spectrogram_peakamounts{i},[0,1,2,3,4,5,6,7]); xlim([1 8]); title({['syllable ' make_temp.parameters.syllables(i) '. Histogram of spectrogram peak amounts ("7" is really 7+)'], ['threshold is' num2str(peak_statistics.peak_height)]}); ylabel('amount of syllables');
% %         saveas(figure(gcf),[save_dir '/peaks_histogram_' make_temp.parameters.syllables(i)], 'fig') 
     
    
end


