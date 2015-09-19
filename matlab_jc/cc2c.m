function vals=cc2c(Predict,Prediction)
for i=[1:14 18:22]
    i
    [a,c]=max(Prediction(i).data(Predict(i).onset:Predict(i).offset));
    b=c+Predict(i).onset;
    for j=1:500
        if (b-j)>Predict(i).onset && (b+j)<Predict(i).offset
            vals(i).data(j)=(Prediction(i).data(b-j)+Prediction(i).data(b+j))/2;
        else
            if (b-j)>Predict(i).onset
                vals(i).data(j)=Prediction(i).data(b-j);
            else
            if (b+j)<Predict(i).offset
                vals(i).data(j)=Prediction(i).data(b+j);
            end
            end
        end
    end
end
