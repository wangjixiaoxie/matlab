% Generate a template and use synshift if you want to

%% Choosing batch file and syllables; making variables

make_temp.parameters.batch_name = input('What is the name of the batch file?  ', 's');
make_temp.parameters.birdname = input('What is the name of the bird? (bk90bk46)  ', 's');
make_temp.parameters.syllables = input('What are the syllables? (no spaces)  ', 's');

make_temp.spec = cell(1,length(make_temp.parameters.syllables));
make_temp.time_range = cell(1,length(make_temp.parameters.syllables));
make_temp.freq_range = cell(1, length(make_temp.parameters.syllables));
make_temp.template = cell(1, length(make_temp.parameters.syllables));
make_temp.parameters.slices = cell(1,length(make_temp.parameters.syllables));

%% Get AVN

% Makes two spectrograms around your syllable of interest. Used to choose
% the points for making your slice (template)
for i = 1:length(make_temp.parameters.syllables)
    [make_temp.spec{i},make_temp.time_range{i},make_temp.freq_range{i}] = get_avn(make_temp.parameters.batch_name,...
        make_temp.parameters.syllables(i),...
        0.2,0.2,'','','obs0');
    
    figure
    imagesc(make_temp.time_range{i}*1000,make_temp.freq_range{i},log(make_temp.spec{i}));syn;ylim([0,1e4]);
    title(['Syllable ' make_temp.parameters.syllables(i)])
    ylabel('Frequency (Hz)')
    xlabel('Time (ms)')
    
    figure
    imagesc(log(make_temp.spec{i}));syn;
    title(['Syllable ' make_temp.parameters.syllables(i)])
    ylabel('Frequency (Hz)')
    xlabel('Points')    
end

%% Choose slices to pull template from and make template

for i = 1:length(make_temp.parameters.syllables)
    % asks for slice points to construct template
    make_temp.parameters.slices{i} = input(['What are the slice coordinates for ' ...
        make_temp.parameters.syllables(i) '?\n(i.e. 60 or 70:82)  ']);
    
    % makes template
    make_temp.template{i} = mean(make_temp.spec{i}(1:2:256,make_temp.parameters.slices{i}),2);
    size(make_temp.template{i}) %should be 128 point vector
    make_temp.template{i}(1:6) = 0;
    make_temp.template{i} = make_temp.template{i}./max(make_temp.template{i});
    mk_tempf(make_temp.parameters.batch_name,make_temp.template{i},2,'obs0');
    
    % figure to look at template
    figure, plot(make_temp.template{i}), title(['syllable ' make_temp.parameters.syllables(i)])
    
    % asks if you want to use syn_shift to add templates predicting how the
    % syllable will shift in frequency
    make_temp.parameters.syn_desire = input(['Would you like to use syn shift on ' make_temp.parameters.syllables(i) '? \n(y or n)  '], 's');

    if strcmpi(make_temp.parameters.syn_desire,'y') == 1
        make_temp.parameters.shift = input(['How much of a shift for ' make_temp.parameters.syllables(i) '?\n(integers: -2, -1, 1:1:4)  ']);
        make_temp.parameters.range = input(['What is the range of the first harmonic of ' make_temp.parameters.syllables(i) '? \n[min max]  ']);
        make_temp.parameters.npeaks = input(['How many peaks are there for ' make_temp.parameters.syllables(i) '?  ']);
        for jj = 1:length(make_temp.parameters.shift)
            make_temp.outvec=zeros(length(make_temp.template{i}),1);
            for ii=1:make_temp.parameters.npeaks
                make_temp.outvecind= ii*(make_temp.parameters.range+make_temp.parameters.shift(jj));
                make_temp.invecind= ii*make_temp.parameters.range;
                make_temp.outvec(make_temp.outvecind(1):make_temp.outvecind(2))=...
                    make_temp.template{i}(make_temp.invecind(1):make_temp.invecind(2),1);
            end
            make_temp.template{i}(:,jj+1) = make_temp.outvec;
        end
    elseif strcmpi(make_temp.parameters.syn_desire, 'n') == 1
    end
    
    % figure to look at template plus synshifts
    figure, plot(make_temp.template{i}), title(['syllable ' make_temp.parameters.syllables(i)])
    
    % writes and saves template to .dat file
    wrt_templ([make_temp.parameters.birdname make_temp.parameters.syllables(i) '.dat'],make_temp.template{i});
end


%% To check use uievtafsim

% for i = 1:length(syllables)
%     try
%         uievtafsim(batch_name,syllables(i))
%     catch err
%         continue
%     end
% end
