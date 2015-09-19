fs=44100;
load '/cobain4/twarren/g100o55/stim/stim3.mat' corpshiftednormg xlist

syllswitch=[14 35]
wavstring2={'9024syll1-7-30.wav' '9024syll2-7-30.wav' '9024syll3-7-30.wav' '9024syll4-7-30.wav' '9024syll5-7-30.wav' '9024syll6-7-30.wav'  }

syll{1}=[xlist(syllswitch(1)) xlist(syllswitch(1)+1)];
syll{2}=[xlist(syllswitch(2)) xlist(syllswitch(2)+1)];
%syll{3}=[xlist(50) xlist(51)];
  maskwid=.005
        
            datnorm_use=corpshiftednormg{1};
            datup_use=corpshiftednormg{3};
            datdown_use=corpshiftednormg{2};
            switchsong{1}=datdown_use;
            switchsong{2}=datup_use;
            switchsong{3}=datnorm_use;
            for i=1:2
                mask{i}=mk_mask(datnorm_use,fs,syll{i}(1),syll{i}(2),maskwid);
            end
count=1;
  for firstswitch=1:3
    for secondswitch=1:2
       
            
            
            
            
           
    
    
    
    
    
    %dat_shift is the mask, to keep previous pitch changes
            dat_shift=datnorm_use.*(1.0-mask{1}) + switchsong{firstswitch}.*mask{1};
            dat_shift=dat_shift.*(1.0-mask{2}) + switchsong{secondswitch}.*mask{2};
            %dat_shift=dat_shift.*(1.0-mask{3}) + switchsong{thirdswitch}.*mask{3};
            datsyll{count}=dat_shift;
            
            %[sm,sp,t,f]=evsmooth(dat_upsngl,fs,0.01);
    %[smnormu,spnormu,tnormu,fnormu]=evsmooth(datnorm_use,fs,0.01);
    %figure(3);clf;ax=zeros([2,1]);
    %ax(1)=subplot(211);imagesc(tnormu,fnormu,log(abs(spnormu)));syn;ylim([0,1e4]);
    %grid on;
    %ax(2)=subplot2(212);imagesc(t,f,log(abs(sp)));syn;ylim([0,1e4]);
    %grid on;
    %linkaxes(ax);
          count=count+1;
        end
    end
  
  
  %for ii=1:2
   %   datsyll{i+8}=datnorm_use.*(1-mask{2}) +switchsong{ii}.*mask{2};
    %  datsyll{i+8}=datsyll{i+8}.*(1-mask{3})+switchsong{ii}.*mask{3}; 
  %end
  for i=1:6
    outwav=strcat('syl',num2str(i),'_',num2str(syllswitch(1)),'_',num2str(syllswitch(2)),'.wav')
    %resample(datsyll{i},44100,44150);
    wavwrite(datsyll{i},44100,16,outwav);
end
  