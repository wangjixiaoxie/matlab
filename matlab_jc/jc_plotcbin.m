
clear fnvl
clear bnds
a='9397_391964989_4_24_11_58_28.wav'; % Based on the first column of the onset matrix (i.e. the output of jc_findonsets).
start1=5.9;  % Based on the second column of the onset matrix (i.e. the output of jc_findonsets).
length1=0.5; % Set this once at the beginning.  Make it the same for all instances of the note.  
end1=start1+length1;
fnvl{1}=a;
bnds{1}=[start1 end1];

for ii=1:length(fnvl)
   %subplot(length(fnvl),1,ii)

   [plainsong,fs]=wavread(a);

   plainsong=floor(plainsong(bnds{ii}(1)*fs:bnds{ii}(2)*fs));
   [sm,sp,t,f]=evsmooth(plainsong,fs,1); % the last variable (after fs) adjusts the amplitude/power threshold for plotting.  Higher values mean higher thresholds.
 
   colormap('jet')
   box off;
   axis on;   
end
