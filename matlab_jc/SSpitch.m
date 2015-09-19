function [peak_pinterp]=SSpitch(fvals,bird_name,syllable)
Fs=32000;
[f_cutoff,t_assay,spect_params]=syllable_params_by_bird(bird_name,syllable);
for ct=1:length(fvals)   
        syl_wav=fvals(ct).datt;
        [S1,F1,T1,P1] =spect_from_waveform(syl_wav,Fs,0,spect_params);
        %disp(['Time to run spect - ' num2str(toc)])
        P1_uncut=P1;
        F1_uncut=F1;
        [r,c]=size(P1);
        whole_P1(1:r,1:c,ct)=P1;
        whole_F1{ct}=F1;
%         if length(T1)>length(save_T1)
%             save_T1=T1;
%         end
        % TO BE REPLACED BY FUN
        f_cut_id=find(F1>f_cutoff(1) & F1<f_cutoff(2));
        F1=F1(f_cut_id);
        P1=P1(f_cut_id,:);


        t_id=find(abs((T1-t_assay))==min(abs(T1-t_assay)));
        spect_slice=P1(:,t_id);
        if size(spect_slice,2)>1
            spect_slice=mean(spect_slice');
        end
        slice_save{ct,1}=spect_slice;

        %figure(2);clf;plot(F1_uncut,P1_uncut(:,t_id));return
        % a measure of amplitude - sum power spectrum
        if length(t_id)>1;spc=mean(sum(P1_uncut(:,t_id)));
        else spc=sum(P1_uncut(:,t_id));end
        sum_power(ct,1)=spc;
%        plot(whole_F1{ct},P1_uncut(:,t_id)),return

        F1_save{ct}=F1;
        max_p_id=find(spect_slice==max(spect_slice));
        peak(ct,1)=F1(max_p_id);
        

        if sum(max_p_id==[1 length(spect_slice)])
        if max_p_id==1
            disp('Peak at LOWER extreme of range')
        elseif max_p_id==length(spect_slice)
            disp('Peak at UPPER extreme of range')
        end
            peak_pinterp(ct,1)=peak(ct,1);
        else
            x_pinterp=F1(max_p_id-1:max_p_id+1);
            y_pinterp=P1(max_p_id-1:max_p_id+1,t_id);
            [peak_pinterp(ct,1),tmp] = pinterp(x_pinterp, y_pinterp);
            if peak_pinterp(ct,1)==0
                disp(['zero here     ' fn])
            end
        end
end
