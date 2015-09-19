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
make_temp.parameters.syllables = input('What are the syllables? (no spaces)  ', 's');

make_temp.spec = cell(1,length(make_temp.parameters.syllables));
make_temp.time_range = cell(1,length(make_temp.parameters.syllables));
make_temp.freq_range = cell(1, length(make_temp.parameters.syllables));
make_temp.template = cell(1, length(make_temp.parameters.syllables));
make_temp.parameters.slices = cell(1,length(make_temp.parameters.syllables));

% dir for saving
suffix=fix(clock);
mkdir(['lt_count_harmonics_results' num2str(suffix(1)) num2str(suffix(2)) num2str(suffix(3)) '_' num2str(suffix(4)) num2str(suffix(5))]);
save_dir=[pwd '/lt_count_harmonics_results' num2str(suffix(1)) num2str(suffix(2)) num2str(suffix(3)) '_' num2str(suffix(4)) num2str(suffix(5))];
 

%% Added by LT, in order to look specifically at syallables followed or preceded by other specific syllables.  (e.g. only look at b if preceded by a)

make_temp.parameters.prent=input('Enter syllables that must come before (blank if none): ', 's');
make_temp.parameters.postnt=input('Enter syllables that must come after (blank if none): ', 's');

%% Get AVN

% Makes two spectrograms around your syllable of interest. Used to choose
% the points for making your slice (template)
for i = 1:length(make_temp.parameters.syllables);
    [make_temp.spec{i},make_temp.time_range{i},make_temp.freq_range{i}] = get_avn(make_temp.parameters.batch_name,...
        make_temp.parameters.syllables(i),...
        0.2,0.2,make_temp.parameters.prent,make_temp.parameters.postnt,'obs0');
    
    figure
    imagesc(make_temp.time_range{i}*1000,make_temp.freq_range{i},log(make_temp.spec{i}));syn;ylim([0,1e4]);
    title(['Syllable ' make_temp.parameters.syllables(i) ', preceded by: ' make_temp.parameters.prent '; followed by: ' make_temp.parameters.postnt ])
    ylabel('Frequency (Hz)')
    xlabel('Time (ms)')
    
    figure
    imagesc(log(make_temp.spec{i}));syn;
    title(['Syllable ' make_temp.parameters.syllables(i) ', preceded by: ' make_temp.parameters.prent '; followed by: ' make_temp.parameters.postnt ])
    ylabel('Frequency (Hz)')
    xlabel('Points')    
end

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

    
    %plot all arbitrary temporary slice
        figure, plot(make_temp.all_spectrogram_slices{i}), title(['syllable ' make_temp.parameters.syllables(i) ' a single instance.'])
%         saveas(figure(gcf), [pwd '/all_spectrogram_slices_' make_temp.parameters.syllables(i)],'fig')
        
        %plot 1 out of every 10 slices
        figure, plot(make_temp.all_spectrogram_slices{i}(:,1:10:size(make_temp.all_spectrogram_slices{i},2))), title(['syllable ' make_temp.parameters.syllables(i) ' a single instance.'])
%         saveas(figure(gcf), [pwd '/all_spectrogram_slices_' make_temp.parameters.syllables(i) '_skipping10'],'fig')

        % plot 1 out of every 25 slices
        figure, plot(make_temp.all_spectrogram_slices{i}(:,1:25:size(make_temp.all_spectrogram_slices{i},2))), title(['syllable ' make_temp.parameters.syllables(i) ' a single instance.'])
%         saveas(figure(gcf), [pwd '/all_spectrogram_slices_' make_temp.parameters.syllables(i) '_skipping25'],'fig')

        
        % find peaks in slices   
        peak_statistics.peak_height=input('what is minimum peak height (range 0-1)? ');   
    for ind=1:size(make_temp.all_spectrograms{i},3);
    
    [peaks,locations]=findpeaks(make_temp.all_spectrogram_slices{i}(:,ind),'MINPEAKHEIGHT',peak_statistics.peak_height);
        make_temp.all_spectrogram_peakheights{i}{ind}=peaks;
        make_temp.all_spectrogram_peaklocations{i}{ind}=locations;
        make_temp.all_spectrogram_peakamounts{i}(ind)=length(peaks);
        
      
    end
            % getting rel amplitudes of 1st 3 harmonics
    %define harmonics (based on inspection of all slices overlayed in an open figure)
            harm=input('what x-coordinates bound your harmonics? (enter "d" for default: [22 32; 47 59; 76 84])');
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
        % plot mean and std for each harmonic (across syllables)
        figure; errorbar([1 2 3], make_temp.all_spectrogram_mean_ampl_three_harmonics{i},make_temp.all_spectrogram_STD_ampl_three_harmonics{i});
        
%         saveas(figure(gcf),[save_dir '/harmonics_mean_and_std', 'fig']) 
       
    
    %plot histogram of peak amounts
        figure;hist(make_temp.all_spectrogram_peakamounts{i},[0,1,2,3,4,5,6,7]); xlim([1 8]); xlabel('number of spectrogram peaks ("7" is really 7+)'); ylabel('amount of syllables');
%         saveas(figure(gcf),[save_dir '/peaks_histogram', 'fig']) 
        disp('press anything to continue')
        
       
     
        
end

% following script: for each syllable, calculate spectral mean using the three harmonics (i.e. freq that is
% midpoint of the three harmonics, weighted by their amplitudes.)
% (e.g. gives range from 1-3 depending on which harmonic dominates.

if(0)
lt_calculate_harmonic_spectral_mean;
end

% save([save_dir '/peak_statistics_' num2str(suffix(1)) num2str(suffix(2)) num2str(suffix(3)) '_' num2str(suffix(4)) num2str(suffix(5))] ,'peak_statistics')
% save([save_dir '/peak_data_' num2str(suffix(1)) num2str(suffix(2)) num2str(suffix(3)) '_' num2str(suffix(4)) num2str(suffix(5))], 'make_temp');

make_template_yesorno = input('do you want to sort through individual slices and make a template? (y or n (=quit and save)) ','s');


%% make a template for syllables passing certain threshold of peak number.


if make_template_yesorno=='y';

for i = 1:length(make_temp.parameters.syllables);
    for j=1:5;
    % getting specs and slices sorted by amount of peaks.
    
    make_temp.syllables_sorted_by_peak_amount{j}.all_spectrograms{i}=make_temp.all_spectrograms{i}(:,:,(make_temp.all_spectrogram_peakamounts{i}==j));
    make_temp.syllables_sorted_by_peak_amount{j}.all_slices{i}=make_temp.all_spectrogram_slices{i}(:, (make_temp.all_spectrogram_peakamounts{i}==j));
    make_temp.syllables_sorted_by_peak_amount{j}.all_peak_locations{i}=make_temp.all_spectrogram_peaklocations{i}(make_temp.all_spectrogram_peakamounts{i}==j);  
    end
    for j=6;
            make_temp.syllables_sorted_by_peak_amount{j}.all_spectrograms{i}=make_temp.all_spectrograms{i}(:,:,(make_temp.all_spectrogram_peakamounts{i}>=j));
    make_temp.syllables_sorted_by_peak_amount{j}.all_slices{i}=make_temp.all_spectrogram_slices{i}(:, (make_temp.all_spectrogram_peakamounts{i}>=j));
    make_temp.syllables_sorted_by_peak_amount{j}.all_peak_locations{i}=make_temp.all_spectrogram_peaklocations{i}(make_temp.all_spectrogram_peakamounts{i}>=j);  

     
    end
end
clear all_peak_locations_compiled
for i=1:length(make_temp.parameters.syllables);
    for ii=1:6
    %plot slices
        figure, plot(make_temp.syllables_sorted_by_peak_amount{ii}.all_slices{i}); title(['All slices with ' num2str(ii) ' peaks.'])

        make_temp.syllables_sorted_by_peak_amount{ii}.mean_slice{i}=squeeze(mean(make_temp.syllables_sorted_by_peak_amount{ii}.all_slices{i},2));
        figure; plot(make_temp.syllables_sorted_by_peak_amount{ii}.mean_slice{i}); title(['mean slice, for notes with ' num2str(ii) ' peaks'])
    
    make_temp.syllables_sorted_by_peak_amount{ii}.mean_spectrogram{i}=squeeze(mean(make_temp.syllables_sorted_by_peak_amount{ii}.all_spectrograms{i},3));
    figure
    imagesc(make_temp.time_range{i}*1000,make_temp.freq_range{i},log(make_temp.syllables_sorted_by_peak_amount{ii}.mean_spectrogram{i}));syn;ylim([0,1e4]);
    title(['Syllable ' make_temp.parameters.syllables(i) ', average for notes with ' num2str(ii) ' peaks.'])
    ylabel('Frequency (Hz)')
    xlabel('Time (ms)')
    
    % plot histogram of average peak locations
    
    if ii~=6;
    all_peak_locations_compiled{i}{ii}=cell2mat(make_temp.syllables_sorted_by_peak_amount{ii}.all_peak_locations{i});
        
    temp=reshape(all_peak_locations_compiled{i}{ii},1,numel(all_peak_locations_compiled{i}{ii}));
    hist(temp,1:2:128); title(['all peak locations combiled from peak class ' num2str(ii)])
    clear temp
    end
    
    %figure showing individual slices.  will show 10, then pause, and
    %wait for command - allow you to note down which ones are outliers and
    %exclude those from the template you will make.
    end
    for ii=1:6;
    for jj=1:size(make_temp.syllables_sorted_by_peak_amount{ii}.all_slices{i},2);
        plot(make_temp.syllables_sorted_by_peak_amount{ii}.all_slices{i}(:,jj)); title(['class of ' num2str(ii) ' peaks.' ' note # ' num2str(jj) 'out of ' num2str(size(make_temp.syllables_sorted_by_peak_amount{ii}.all_slices{i},2))])
        instruction=input('press ENTER to go to next note. press "2" to skip 10 notes. press "3" to skip to next peak class');
        if instruction ==2;
            jj=jj+10;
        end
        if instruction ==3;
            
        break        
        end
    end
    
    end
    

%make the template .dat file. 

% % if(0) input('any notes from classes 1-4 that you wish to include in template (y or n)?', 's')=='y';
% %     for ii=1:6;
% %     omitted_slices{ii}=input
% %     desired_slices{ii}=input('which
% %     end
% % end

    clear temp
if(0)
    
    for ii=5;
    omitted_slices{ii}=input('to make a template from class 4,5, and 6 combined, which slices do you want to remove (e.g.[2 3 24]?');
    omitted_slices{ii}
%     while input('do omitted slices look correct (y or n)')=='n';
%             omitted_slices{ii}=input('to make a template from class 5, which slices do you want to remove (e.g.[2 3 24]?');
%     end
    temp.slices = make_temp.syllables_sorted_by_peak_amount{ii}.all_slices{i};
    temp.spectrograms=make_temp.syllables_sorted_by_peak_amount{ii}.all_spectrograms{i};
    temp.slices(:,omitted_slices{ii})=[];
    
    temp.spectrograms(:,:,omitted_slices{ii})=[];
    
    size(temp.slices,2), size(make_temp.syllables_sorted_by_peak_amount{ii}.all_slices{i},2)
   
    temp.mean_spectrogram=mean(temp.spectrograms,3);
    imagesc(make_temp.time_range{i}*1000,make_temp.freq_range{i},log(temp.mean_spectrogram));syn;ylim([0,1e4]);
    title(['final template'])
    ylabel('Frequency (Hz)')
    xlabel('Time (ms)')

    
    temp.templ=mean(temp.slices,2);
    temp.templ=temp.templ./max(temp.templ);
    
    
    %plot the final template
    plot(temp.templ); title('final template')
    
    %making the template .dat file
    wrt_templ([make_temp.parameters.birdname make_temp.parameters.syllables(i) '_' num2str(ii) 'peaks_' num2str(peak_height) 'tshld' '.dat'],temp.templ);


end

end
end

%pick out the slices that you want by hand

clear desired_slices
clear temp
clear dummy


desired_slices{3}=input('enter the notes you want in class 3 (e.g. [3:6 7 8 10:16]');
desired_slices{4}=input('enter the notes you want in class 4 (e.g. [3:6 7 8 10:16]');
desired_slices{5}=input('enter the notes you want in class 5 (e.g. [3:6 7 8 10:16]');
desired_slices{6}=input('enter the notes you want in class 6 (e.g. [3:6 7 8 10:16]');


for ii=3:6;
temp.slices{ii} = make_temp.syllables_sorted_by_peak_amount{ii}.all_slices{i};
temp.spectrograms{ii}=make_temp.syllables_sorted_by_peak_amount{ii}.all_spectrograms{i};
temp.slices{ii}=temp.slices{ii}(:,desired_slices{ii});
temp.spectrograms{ii}=temp.spectrograms{ii}(:,:,desired_slices{ii});


  
    size(temp.slices{ii},2), size(desired_slices{ii})
    disp('verify that the quantity of desired slices and selected slices are equal')
    pause
   
end

temp.all_desired_slices=cell2mat(temp.slices);

for ii=3:6
    if ii==3;
        temp.all_desired_specs=temp.spectrograms{ii};
        continue
        
    end
    if ii~=3;
        temp.all_desired_specs=cat(3, temp.all_desired_specs, temp.spectrograms{ii});
    end
end

            



temp.mean_all_spectrogram=mean(temp.all_desired_specs,3);

    
    %plot
    figure; imagesc(make_temp.time_range{i}*1000,make_temp.freq_range{i},log(temp.mean_all_spectrogram));syn;ylim([0,1e4]);
    title(['final template spectrogram'])
    ylabel('Frequency (Hz)')
    xlabel('Time (ms)')
pause
    
    temp.templ=mean(temp.all_desired_slices,2);
    temp.templ=temp.templ./max(temp.templ);
    
    
    %plot the final template
    figure; plot(temp.templ); title('final template')
    
    %making the template .dat file
    wrt_templ([make_temp.parameters.birdname make_temp.parameters.syllables(i) '_classes3to6combined_' num2str(peak_height) 'tshld' '.dat'],temp.templ);

 
    
%saving things


save([save_dir '/peak_statistics_' num2str(suffix(1)) num2str(suffix(2)) num2str(suffix(3)) '_' num2str(suffix(4)) num2str(suffix(5))] ,'peak_statistics')
save([save_dir '/peak_data_' num2str(suffix(1)) num2str(suffix(2)) num2str(suffix(3)) '_' num2str(suffix(4)) num2str(suffix(5))], 'make_temp');

end
