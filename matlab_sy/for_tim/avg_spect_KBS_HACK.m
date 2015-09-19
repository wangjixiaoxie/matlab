function [avg_spect, std_spect, axisdata]=avg_spect(batch_in);

%return average and std of spectrograms read in from batch file of idx spects
%if plot flag ==1 make plot of matrices
%
% if norm_flag == 1 then normalize each spectrogram
%    want normalization to be robust to 1) high amplitude outliers
%                                   and 2) differing amounts of background (vs song)
%    so: throw out background (for reason #2) set by norm_floor
%        and take median value for normalization (less senstive to outliers)
%axisdata= [t_min, t_max, f_min, f_max];


norm_flag=0;
norm_floor=85;  %number of decibels below max that is excluded from spect before norm
plot_flag=1;


%open batch_file
meta_fid=fopen([batch_in]);
if meta_fid == -1 | batch_in == 0
      disp('cannot open file' )
      disp (batch_in)
      return
end

n_in=0;  %to accumulate number of spectrograms

while 1 
     %get filename
     spectfile = fscanf(meta_fid,'%s',1);
     %end when there are no more spectfiles 
     if isempty(spectfile);
        break
     end
     
      %if spectfile exists, get it
     if exist([spectfile]);   
       [spect, nfft, spect_win, noverlap, t_min, t_max, f_min, f_max] = read_spect(spectfile);
       [n_rows,n_cols]=size(spect);
       %flip spectrogram so that 100db is loudest, 0db softest
       spect=100-spect;
       n_in=n_in+1
     else
       disp(['cannot find ', spectfile])
     end
     
     
     % normalize spectrograms here if desired
     if norm_flag==1 
       % now collect values that are above the norm_floor and get their median
       spect_vect=reshape(spect,[n_rows*n_cols, 1]);    
       keep_ind=find(spect_vect>=(100-norm_floor));
       norm_val=median(spect_vect(keep_ind));     
       %subtract median value from spect so that all values are now relative to median
       spect=spect-norm_val;
     end
     
     %accumulate spectrograms in matrix
     
     if n_in==1;
       all_spect=[spect];
       [n_rows_all,n_cols_all]=size(spect);
       max_cols=n_cols_all;
       min_cols=n_cols_all;
     else
        if n_rows~=n_rows_all;
	  disp('different number of rows...abort!') 
	  break 
	end   
	if n_cols~=n_cols_all 
	  if n_cols > n_cols_all
	    spect=spect(:,1:n_cols_all);
	    max_cols=max(max_cols,n_cols);	  
	  else
            all_spect=all_spect(:,1:n_cols);	    
	    n_cols_all=n_cols;
	    min_cols=min(min_cols,n_cols);
	  end	   
        end
        all_spect=all_spect+spect;
     end

end 


avg_spect=all_spect./n_in;

%now get std dev
%frewind here and accumulate average squared deviations from mean to get std 
std_spect=[];
frewind(meta_fid);

n_in=0;
[n_rows_all,n_cols_all]=size(all_spect);
while 1 
     %get filename
     spectfile = fscanf(meta_fid,'%s',1);
     %end when there are no more spectfiles 
     if isempty(spectfile);
        break
     end
     
      %if spectfile exists, get it
     if exist([spectfile]);   
       [spect, nfft, spect_win, noverlap, t_min, t_max, f_min, f_max] = read_spect(spectfile);      
       [n_rows,n_cols]=size(spect);
       %flip spectrogram so that 100db is loudest, 0db softest
       spect=100-spect;
       n_in=n_in+1
     else
       disp(['cannot find ', spectfile])
     end
     
     
     % normalize spectrograms here if desired
     
     if norm_flag==1        
       % now collect values that are above the norm_floor and get their median
       spect_vect=reshape(spect,[n_rows*n_cols, 1]);    
       keep_ind=find(spect_vect>=(100-norm_floor));
       norm_val=median(spect_vect(keep_ind));     
       %subtract median value from spect so that all values are now relative to median
       spect=spect-norm_val;
     end
     
    
     %take difference from mean; now just cut to mean spect size
         
     [n_rows,n_cols]=size(spect);
     if n_rows~=n_rows_all;
	disp('different number of rows...abort!') 
	break 
     end   
     if n_cols~=n_cols_all 
	spect=spect(:,1:n_cols_all);
     end
     diff_spect=spect-avg_spect;
     diff_spect_sqr=diff_spect.^2;  %square deviations
     if n_in==1
       std_spect=diff_spect_sqr;
     else
       std_spect=std_spect+diff_spect_sqr;
     end  
end 
%divide by n-1
std_spect=std_spect./(n_in-1);
%take sqrt
std_spect=std_spect.^.5;

fclose(meta_fid);
disp(['columns size ranged from ',num2str(min_cols),' to ',num2str(max_cols)]);
disp(['number of spectrograms = ',num2str(n_in)]);

axisdata=[t_min,t_max,f_min,f_max];

%flip spectrogram frequency axis
avg_spect=flipud(avg_spect);
std_spect=flipud(std_spect);

if plot_flag==1
  figure
  subplot(2,1,1)
  imagesc(avg_spect);
  ylabel('mean');
  subplot(2,1,2)
  imagesc(std_spect);
  ylabel('standard deviation');
  subplot(2,1,1)
end
