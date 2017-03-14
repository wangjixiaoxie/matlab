function [histbins allcounts] = mbatchampdist(batch)
%
% returns distribution of sound amplitude values accumulated across all
% sound files in batch.
%
%
%
%

fid=fopen(batch,'r');

histbins = [-12:.05:12]; % hardcoded freq. bins for histogram.
allcounts = zeros(1,length(histbins));
%allvals = zeros(1,500000);
i=1;
%disp('mBatchHistAmp() working...');
while (1)
    fn=fgetl(fid);
    if (~ischar(fn))
        break;
    end
    if (~exist(fn,'file'))
        continue;
    end
    
    if (~mod(i,20)) %every 20th file
        
        %disp(fn);
        
        [pth,nm,ext]=fileparts(fn);
        if (strcmp(ext,'.ebin'))
            [dat,fs]=readevtaf(fn,'0r');
            sm=evsmooth(dat,fs,0.01);
        elseif(strcmp(ext,'.cbin'))
            [dat,fs]=ReadCbinFile(fn);
            if size(dat,2)>1 % extract audio channel
                dat = dat(:,1);
            end
            sm=mquicksmooth(dat,fs);
        elseif(strcmp(ext,'.wav'))
            [dat,fs]=wavread(fn);
            sm=mquicksmooth(dat,fs);
        elseif(strcmp(ext, '.rhd'))
            [frequency_parameters, board_adc_data] = pj_readIntanNoGui_AudioOnly(fn);
            fs=frequency_parameters.amplifier_sample_rate; % for neural, analog, and dig
            dat = board_adc_data(1,:);
            sm=mquicksmooth(dat,fs);
        end
        
        sm = log10(sm);
        
        %    if(sum(allvals)==0)
        %
        %        %allvals = log10(sm);
        %    else
        %        %allvals = cat(1,allvals,log10(sm));
        %    end
        [theHist] = hist(sm,histbins);
        allcounts = allcounts+theHist;
    end
    i=i+1;
end

fclose(fid);
%allcounts = hist(allvals,histbins);
%disp('mBatchHistAmp() done');

return;