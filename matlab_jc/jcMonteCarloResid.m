function [zSCORE]=jcMonteCarloResid(normA1,normB1,normA2,normB2,a,b,window_width)
the_length=size(normA1,2);

normA1=normA1';
normB1=normB1';
normA2=normA2';
normB2=normB2';
    
    
%UD_PRE
count1=1;
count2=1;
for i=1:the_length
    m(1,i)=median(normA1(i,a:a+window_width));
    m(2,i)=median(normB1(i,b:b+window_width));
    count1=count1+1;
end

%D_PRE
for j=1:the_length
    m(3,j)=median(normA2(j,a:a+window_width));
    m(4,j)=median(normB2(j,b:b+window_width));
    count2=count2+1;
end

ccLMAN=corrcoef(m(1,:),m(2,:));
ccLMAN=ccLMAN(2);
ccnone=corrcoef(m(3,:),m(4,:));
ccnone=ccnone(2);
for k=1:4
    stand(k)=std(m(k,:));
end

%To avoid imaginary answers
if stand(1)>stand(3)
    sdLMAN1=sqrt(stand(1)^2-stand(3)^2);
else sdLMAN1=0;
end
if stand(2)>stand(4)
    sdLMAN2=sqrt(stand(2)^2-stand(4)^2);
else sdLMAN2=0;
end

if stand(4)>stand(2) && stand(3)>stand(1)
    zSCORE=100;
else


%The Monte Carlo simulation
for ii=1:1000
    noisy1=randn(1,the_length)*sdLMAN1;
    noisy2=randn(1,the_length)*sdLMAN2;
    mfake1=noisy1+m(3,:); %mocks m(1,:)
    mfake2=noisy2+m(4,:); %mocks m(2,:)
    cc=corrcoef(mfake1,mfake2);
    ccMock(ii)=cc(2);
end

zSCORE=(ccLMAN-mean(ccMock))/(std(ccMock));
end