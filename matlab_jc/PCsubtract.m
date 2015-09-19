function Reconst=PCsubtract(condition,SCORE,PC,LATENT,ACsize,INAsize)
if condition==1
% Simulate LMAN only
    for i=1:100
            clear PCmock
            % For each principal component
            for k=1:size(SCORE,2)
                % pick random ACSF pitch curve
                g=round(rand*(ACsize-1)+1);
                % pick random INA pitch curve
                h=round(rand*(INAsize-1)+1);
                % Subtract difference in magnitude of score/weight of PC #k for
                % the two randomly selected curves
                aa=(abs(SCORE(g,k))-abs(SCORE(h+(ACsize),k)));
                % Randomize sign of the simulated PC difference
                if rand>0.5
                    sign=1;
                else
                    sign=-1;
                end
                PCmock(:,k)=abs(aa)*PC(:,k);
            end
            Reconst(:,i)=sum(PCmock');
    end
else 
    if condition==2
        % Simulate Both only
        for i=1:100
                clear PCmock
                % For each principal component
                for k=1:size(SCORE,2)
                    % pick random ACSF pitch curve
                    g=round(rand*(ACsize-1)+1);
                    % pick random INA pitch curve
                    % Subtract difference in magnitude of score/weight of PC #k for
                    % the two randomly selected curves
                    aa=abs(SCORE(g,k));
                    % Randomize sign of the simulated PC difference
                    if rand>0.5
                        sign=1;
                    else
                        sign=-1;
                    end
                    PCmock(:,k)=abs(aa)*PC(:,k);
                end
                Reconst(:,i)=sum(PCmock');
        end
    else
        if condition==3
                    % Simulate non-LMAN only
        for i=1:100
                clear PCmock
                % For each principal component
                for k=1:size(SCORE,2)
                    % pick random INA pitch curve
                    h=round(rand*(INAsize-1)+1);
                    % Subtract difference in magnitude of score/weight of PC #k for
                    % the two randomly selected curves
                    aa=abs(SCORE(h+(ACsize),k));
                    % Randomize sign of the simulated PC difference
                    if rand>0.5
                        sign=1;
                    else
                        sign=-1;
                    end
                    PCmock(:,k)=abs(aa)*PC(:,k);
                end
                Reconst(:,i)=sum(PCmock');
        end
        end
    end
end

