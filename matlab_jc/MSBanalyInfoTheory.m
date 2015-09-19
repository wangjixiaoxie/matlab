count=0;
for j=0.4:0.1:1
    count=count+1;
    for i=1:28

        if isequal(Predict(i).direction,'up')
            valu=60;
        else
            valu=40;
        end
        if i<23
            [pretime(count,i) posttime(count,i)]=ContingSimMSB(Predict(i).Targeting-Predict(i).onset,Predict(i).ResidAC(Predict(i).onset:Predict(i).onset+400,:),valu,j);
        else
        	[pretime(count,i) posttime(count,i)]=ContingSimMSB(Predict(i).Targeting-Predict(i).onset,Predict(i).ResidAC(Predict(i).onset:Predict(i).offset,:),valu,j);
        end
    end
end

predicted=mean(posttime'-pretime')/8;
for i=1:22
    notewidth(i)=400;
end
for i=23:28
    notewidth(i)=Predict(i).offset-Predict(i).onset;
end
actual=mean(notewidth)/8;

[h,p]=ttest(posttime(5,:)-pretime(5,:),notewidth)
p=2.3e-8