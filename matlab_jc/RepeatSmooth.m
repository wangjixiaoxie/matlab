function [smooth,avgnote]=RepeatSmooth(shifted)
%The raw data array is the output of findwnote2.
fs=32000;
%Smooth the raw data
for i=1:size(shifted,1)
    [holder]=SmoothData(shifted(i,:),fs,1);
    smooth(i,:)=log(holder);
end

%Get average smoothed data
for i=1:size(smooth,2)
    sumnote(i)=0;
    for j=1:size(smooth,1)
        sumnote(i)=sumnote(i)+smooth(j,i);
    end
    avgnote(i)=sumnote(i)/size(smooth,1);
end
% 
% %Align the pitch curves
% for i=1:size(smooth,1)
%     h=xcov(avgnote(1:15000),smooth(i,1:15000));
%     [peak,index]=max(h);
%     shift(i)=(index-15000);
%     starter=500-round(shift(i));
%     ender=starter+14499;
%     aligned(i,:)=[zeros(1,500),shifted(i,starter:ender)];
%     
% end
