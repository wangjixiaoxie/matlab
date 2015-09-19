function [spm,mu_spect,std_spect,axis_data,p_spect,p_snd]=avg_spect_proto(batch_in,syl,proto_sng,p_f_type);
% xcov power vs. time profile of syl with prototypical to center
%
% store spects in matrix
%
% get mean and var of spects point-by-point
%
%
%return average and std of spectrograms read in from batch file of idx spects
%if plot flag ==1 make plot of matrices
%
% if norm_flag == 1 then normalize each spectrogram
%    want normalization to be robust to 1) high amplitude outliers
%                                   and 2) differing amounts of background (vs song)
%    so: throw out background (for reason #2) set by norm_floor
%        and take median value for normalization (less senstive to outliers)
%axisdata= [t_min, t_max, f_min, f_max];

%number of points to pull out pre and post sound ... all sounds are
%resampled at 32e3
char(syl)
pre = 800;
post = 500;


disp('Resamples all sounds at 32k')

norm_flag=0;

% was 90
norm_floor=290;  %number of decibels below max that is excluded from spect before norm
plot_flag=1;

%values for calculating spectrogram
spect_win_dur=10;
spect_overlap = .80;  %percentage of overlap of specgram window
nfft=round(32e3*spect_win_dur/1000);
nfft = 2^nextpow2(nfft);
spect_win = hanning(nfft);
noverlap = round(spect_overlap*length(spect_win)); %number of overlapping points   


% get the prototypical syllable
[p_snd,Fs_p] = soundin('',char(proto_sng),char(p_f_type));
n_syls = 1;
if Fs_p == -1
    Fs_p = 32e3;
    disp('Warning!!! FS assumed to be 32000 for f_type = b')
end

%resample the sound at 32k
p_snd = resample(p_snd,32e3,Fs_p);

%load in batchfile associated with snd file and find syllable
load([char(proto_sng) '.not.mat']);
i = strfind(char(labels),char(syl));
p_snd = p_snd((round(onsets(i)*32)-pre):round(offsets(i)*32)+post)';    
lng_snd = length(p_snd);

%now calculate spectrogram
[spect, freq, time] = specgram(p_snd, nfft, 32e3, spect_win, noverlap);
spect=real(log10(spect));disp('SAM HACK LOG POWER')
spect = spect(find(freq <10000),:); %changed from 8k to match new filtering methods from autolabelkb
freq = freq(find(freq <10000));
%flip spectrogram so that 100db is loudest, 0db softest

%spect=100-spect;  % HACK commented out

%modify spectrogram appropriately ...i.e. normalize 
[n_rows,n_cols]=size(spect);
p_lng = n_cols;
% now collect values that are above the norm_floor and get their median
spect_vect=reshape(spect,[n_rows*n_cols, 1]);    

% HACK COMMENTED OUT
% keep_ind=find(spect_vect>=(100-norm_floor));  %HERE
% norm_val=median(spect_vect(keep_ind));     


%subtract median value from spect so that all values are now relative to median
p_spect=spect;%-norm_val; % HACK COMMENTED OUT

%sum accross frequencies for each time point to get temporal power
%profile
p_spect = flipud(abs(p_spect));
p_spect = resample(p_spect',4,1)';
[p_freq,p_time]=size(p_spect);p_lng = p_time;

% HACK COMMENTED OUT
%inds = p_spect<=0; p_spect(p_spect<=0) = 0;
p_pwr_v_time = sum(p_spect);
autocor_psyl = max(xcorr(p_pwr_v_time));


%values for plotting are based on the prototypical syllable
f_min = freq(1);
f_max = freq(length(freq));
freq_spect = [f_min, f_max];
t_min = time(1)*1000; %convert to ms
t_max = time(length(time))*1000; %convert to ms
%adjust time axis for spectrogram offset (1/2 window duration in ms)
t_min = t_min + 0.5*nfft/32;  
t_max = t_max + 0.5*nfft/32;  
time_spect = [t_min, t_max];    
axis_data = [time_spect freq_spect];

%open batch_file
meta_fid=fopen([batch_in]);
if meta_fid == -1 | batch_in == 0
      disp('cannot open file' )
      disp (batch_in)
      return
end

n_syls = 0;
while 1 
     %get filename
     file = fscanf(meta_fid,'%s',1);
     %end when there are no more spectfiles 

%     if isempty(file) | n_syls>5;
     if isempty(file) | n_syls>100;
         %if isempty(file);
%         if n_syls>5;disp('n_syls>5 - STOPPING EARLY');end
         if n_syls>100;disp('n_syls>100 - STOPPING EARLY');end
         break
     end

     if isempty(strfind(file,'.filt'))
        file = [file '.filt'];
     end
     
     %if file exists, get it
     if exist([file]);   
         [song, Fs] = read_filt(file);
     else
       disp(['cannot find ', file])
     end
     
     if Fs ==-1
         Fs = 32000;
     end
     
     %resample the song
     sng = resample(song,32e3,Fs);
     
     load([char(file(1:end-4)) 'not.mat']);
     
	% get the syl
	for i = 1:length(labels)
        if strcmp(labels(i),char(syl));
            n_syls = n_syls+1;
            if mod(n_syls,100) == 0
                n_syls
            end
            if round(onsets(i)*32)-pre >=1 & post+round(offsets(i)*32) <= length(sng)
                snd = sng(round(onsets(i)*32)-pre:post+round(offsets(i)*32))';
            elseif round(onsets(i)*32)-pre < 1
                snd = sng(1:post+round(offsets(i)*32))';
            elseif post+round(offsets(i)*32) <= length(sng)
                snd = sng(round(onsets(i)*32)-pre:end)';
            else
                snd = [];
            end
            
            if ~isempty(snd)
                %lets time warp the sound with phase-vocoder
                
                %snd = pvoc(snd,round(100*length(snd)/lng_snd)/100,nfft);
                %now calculate spectrogram
                [spect, freq, time] = specgram(snd, nfft, 32e3, spect_win, noverlap);
                spect=real(log10(spect));
                
              spect(find(spect<1))=0;  disp('SAM HACK FLOOR')
              
                spect = spect(find(freq <10000),:);%changed from 8k to match new filtering

                freq = freq(find(freq <10000));
                %modify spectrogram appropriately ...i.e. normalize 
                [n_rows,n_cols]=size(spect);
                %flip spectrogram so that 100db is loudest, 0db softest
%                spect=100-spect; % normalize and subtract to switch
                
                % now collect values that are above the norm_floor and get their median
                spect_vect=reshape(spect,[n_rows*n_cols, 1]);    
                
                % HACK COMMENTED OUT
%                 keep_ind=find(spect_vect>=(100-norm_floor)); %HERE
%                 norm_val=median(spect_vect(keep_ind));     
                %subtract median value from spect so that all values are now relative to median
%                spect=spect-norm_val;   % HACK COMMENTED OUT

                spect = flipud(abs(spect));
                spect = resample(spect',4,1)';
                syl_pwr_v_time = sum(spect);
                szcrnt = size(spect);
                if szcrnt(1) > 1 & szcrnt(2) > 1
                    
                    %lets try doing a linear time warping hear
                    %                 %calculate shift from max of xcorr(p_pwr_v_time,syl_pwr_v_time)
                    rszspect = interp2(spect,1:p_time,[1:p_freq]');
                    %                 %calculate the cross-corralation and auto_correlation
                    xcorvect = xcorr(p_pwr_v_time,syl_pwr_v_time);
                    autocor_syl = max(xcorr(syl_pwr_v_time));
                    %normalize xcorvect by the max of autocorrelations
                    nrmxcorvect = xcorvect./sqrt(autocor_psyl*autocor_syl);
                    %now convert back into FS units by resampling
                    rs_nrm_xcorvect = resample(nrmxcorvect,1,2);
                    %find the max value in a window around center of xcorr vector
                    z_pnt = median([1:length(rs_nrm_xcorvect)]);
                    maxval = max(rs_nrm_xcorvect(round(z_pnt-10):round(z_pnt+10)));
                    pntofmax = find(rs_nrm_xcorvect == maxval);
                    %the lag in points to the center of mass of xcorr
                    lag2cntrM = round((pntofmax-z_pnt));
                    %now shift spect by lag
                    if lag2cntrM >0
                        spect = [spect(:,[lag2cntrM:end]) zeros(n_rows,lag2cntrM)];
                    elseif lag2cntrM <0
                        spect = [zeros(n_rows,(abs(lag2cntrM))) spect(:,[1:end-abs(lag2cntrM)])];
                    end
                    [n_rows,n_cols]=size(spect);
                    %finally, normalize to p_length
                    diff = p_lng-n_cols;
                    if diff > 0
                        pntspre = zeros(n_rows,floor(diff/2));
                        pntspost = zeros(n_rows,diff-floor(diff/2));
                        rszspect = [pntspre spect pntspost];
                    elseif diff < 0
                        rszspect = spect(:,[1+abs(floor(diff/2)):end-(abs(diff)-abs(floor(diff/2)))]);
                    end

                else n_syls = n_syls-1;
                    
                end
                spm(:,:,n_syls) = rszspect;
            end
        end
	end
end

% calculate varSPM
mu_spect = mean(spm,3);

[n_rows,n_cols] = size(mu_spect);
% now collect values that are above the norm_floor and get their median
mu_spect_vect=reshape(mu_spect,[n_rows*n_cols, 1]);    
keep_ind=find(mu_spect_vect>=(100-norm_floor));
norm_val=median(mu_spect_vect(keep_ind));     
%subtract median value from spect so that all values are now relative to median
%mu_spect=mu_spect-norm_val;  % HACK COMMENTED OUT

% HACK COMMENTED OUT - THIS MADE A BIG DIFFERENCE
% KBS- DO TAKE LOG SPECTRUM, BUT FIX THE BELOW TO AVOID NEGATIVE VALS
%
%inds = mu_spect<=0; mu_spect(mu_spect<=0) = 0;

for l = 1:n_syls
    std_spect = sum(spm-repmat(mu_spect,[1 1 n_syls]),3);
end
