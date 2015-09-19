function vals=cc1c(Predict)
for i=1:22
    [a,c]=max(Predict(i).LearnedNorm(Predict(i).onset:Predict(i).offset));
    b=c+Predict(i).onset;
    for j=1:500
        if (b-j)>Predict(i).onset && (b+j)<Predict(i).offset
            vals(i).data(j)=(Predict(i).LearnedNorm(b-j)+Predict(i).LearnedNorm(b+j))/2;
        else
            if (b-j)>Predict(i).onset
                vals(i).data(j)=Predict(i).LearnedNorm(b-j);
            else
            if (b+j)<Predict(i).offset
                vals(i).data(j)=Predict(i).LearnedNorm(b+j);
            end
            end
        end
    end
end
