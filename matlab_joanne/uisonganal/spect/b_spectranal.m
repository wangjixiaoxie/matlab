function [nspect, sn, bf, modul, br, pva, spect_sum, f_centers]=b_spectranal(notefiles, norm_flag, filetype, d_flag);

% need to strip first bin of spectra: see psdanal for model that works
%Note: function deltaprods assumes that first bin of spectrum represents a frequency of f_step, NOT 0.



%set display options
if d_flag == 1
  sum_d_flag=1;
  d_flag=0;
elseif d_flag==2
  sum_d_flag=1;
  d_flag=1;
else
  sum_d_flag=0;
  d_flag=0;
end



 
%initialize variables
sn=[];  %signal to noise
bf=[]; %best fundamental frequency (within range set in deltaprods)- 0 implies maxima at boundary of range
modul=[]; %length of normalized spectrum, see spectranla fro normalization options
br=[]; %best ratio of cos products to sin products
pva=[]; %ratio of amplitude at peaks of spectrum to mean amplitude of spectrum (power instead?)





%open notefiles

   metafile1=notefiles;
   meta_fid1=fopen([notefiles]);
   if meta_fid1 == -1 | metafile1 == 0
     		 disp('cannot open file' )
      		disp (metafile1)
   end
   
   
%until out of files, read a notefile

while 1
  
   %get notefile name
     notefile = fscanf(meta_fid1,'%s',1);
   %end when there are no more notefiles
     if (notefile==[])
        %disp('End of notefiles')
        break
     end
     
     
   %convert notefile name to songfile name
   ind=findstr('.not.mat',notefile);
   if isempty(ind)
      disp(['batchfile has non-notefile entry: ',notefile])
      break
   else
      soundfile=notefile(1:ind-1);
   end
   
   
   %calculate all sorts of stuff
   [nspect1, note_sn, best_freq, modulation, best_ratio, peak_vs_avg, f_centers]=spectranal(soundfile, norm_flag, filetype, d_flag);
   
   %concatenate results
   sn=[sn,note_sn];
   bf=[bf,best_freq];
   modul=[modul,modulation];
   br=[br,best_ratio];
   pva=[pva,peak_vs_avg];
   nspect=[nspect,nspect1];
   
end



spect_sum=[sn;bf;modul;br;pva];


if sum_d_flag == 1
 figure
 plot(pva,bf,'+')
end


   
   
