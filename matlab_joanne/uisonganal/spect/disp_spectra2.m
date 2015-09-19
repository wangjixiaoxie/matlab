function [score,ind1,ind2]=disp_spectra2(spect1,spect2,n, freq_c1, freq_c2,perm_flag);

%2 alternative forced choice comparison of spectra
% first rom of score contains 1 if spect from spect1 was chosen, else 2
% second row contains button that was used = confidence

%default is for perm_flag = 1
if nargin <=5
  perm_flag=1;
elseif perm_flag~=0
  perm_flag=1;
end



%norm flag determines vertical scaling

norm_flag=0;

l1=size(spect1,2);
l2=size(spect2,2);
n=min([n,l1,l2]);

ind1=randperm(l1);
ind2=randperm(l2);

ind1=ind1(1:n);
ind2=ind2(1:n);
score=zeros(2,n);

if perm_flag==1
  top1=round(rand(1,n));
else
  top1=ones(1,n);
end

if norm_flag==1
 max_y1=max([spect1(:);spect2(:)]);
 max_y2=max_y1;
 max_y2=12;
 max_y1=12;
end


figure
for i=1:n 
   disp(['n = ',num2str(i)]); 
  if top1(i)==1
   subplot(2,1,1)
   plot(freq_c1,spect1(:,ind1(i)))   
   if norm_flag==1
     set(gca,'ylim',[0,max_y1]);
   end
   %label only when not permuted
   if perm_flag==0
     title(['top plot from arg1, row ',num2str(ind1(i))]);
   end   
   subplot(2,1,2)  
   plot(freq_c2,spect2(:,ind2(i)))
   if norm_flag==1
     set(gca,'ylim',[0,max_y2]);
   end    
   if perm_flag==0
     title(['from arg2, row ',num2str(ind2(i))]);
   end
   %in=input('top more harmonic = 1; bottom more = 2','s');
   ylim=get(gca,'ylim');
   [x,y,but]=ginput(1); 
   if y > 1.1*ylim(2)
     in = 'top';
   else
     in = 'bottom';
   end  
   if strcmp(in,'top');
     score(1,i)=1;
     disp('top');
   elseif strcmp(in,'bottom');
     score(1,i)=2;
     disp('bottom');
   else
     score(i) == -1;
     disp('funky');
   end
   score(2,i)=but;
  end
  if top1(i)== 0
   subplot(2,1,1)
   plot(freq_c2, spect2(:,ind2(i)))
   if norm_flag==1
     set(gca,'ylim',[0,max_y2]);
   end      
   subplot(2,1,2)
   plot(freq_c1, spect1(:,ind1(i)))
   if norm_flag==1
     set(gca,'ylim',[0,max_y1]);
   end      
   %in=input('top more harmonic = 1; bottom more = 2','s');
   ylim=get(gca,'ylim');
   [x,y,but]=ginput(1); 
   if y > 1.1*ylim(2);
     in = 'top';
   else
     in = 'bottom';
   end     
   if strcmp(in,'top');
     score(1,i)=2;
     disp('top');
   elseif strcmp(in,'bottom');
     score(1,i)=1;
     disp('bottom');
   else
     score(i) == -1;
     disp('funky');
   end
   score(2,i)=but;
  end
end



%display summary

subplot(1,1,1);

%overall
scored1=find(score(1,:)==1);
scored2=find(score(1,:)==2);
nscored1=length(scored1);
nscored2=length(scored2);
per1=100*nscored1/(nscored1+nscored2);
plot([0,4],[per1,per1]);
hold on;

xlabel('confidence');
ylabel(['% arg1 scored less noisy (mean=',num2str(per1),')']);

%by confidence
for i=1:3
 scored1=find(score(1,:)==1 & score(2,:) == i);
 scored2=find(score(1,:)==2 & score(2,:) == i);
 nscored1=length(scored1);
 nscored2=length(scored2);
 per1=100*nscored1/(nscored1+nscored2);
 plot(i,per1,'o');
 text(i-.3,per1+5,[num2str(per1),'% (n=',num2str(nscored1+nscored2),')']);
 axis([0,4,0,100]);
end

if exist('inputname.m')
  title([inputname(1),' vs ',inputname(2)]);
end

subplot(1,1,1);

