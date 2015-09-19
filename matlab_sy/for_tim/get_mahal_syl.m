function get_mahal_syl(templatestrct,batch2srch,name,min)

%function mahalmin = get_mahal_syl(templatestrct,labels2fnd,batch2srch,name)
%
%<templatestrct>: a template structure from autolabel
%<labels>: the labels you want to attach sounds to
%<batch2srch>: the batchfile used in template creation
%<name>:base of the name to save mhlsyl.strct.mat as
%% to add:
%<dist_type>: either max or min;
%max finds a syllable that is distant from the centroid of that type (in the mahalanobus sense)
%min finds the closest syllable to the centroid of a particulr label
% can be thought of as finding an outlier and a prototype syllable,
% respectively

warning('off','all');

%first unravel the template structure

template_mean = templatestrct.template_mean;
template_cov = templatestrct.template_cov;
labelIndex = templatestrct.labelIndex;
for j = 1:length(labelIndex)
    for h = 1:10
        d(j).top10(h) = 1000000;
        mahalmin{j}.syl{h}.snd = [];
        mahalmin{j}.dist(h) = 1000000;
    end
end
%for each song in the file
fidB = fopen(char(batch2srch))
while ~feof(fidB)

    %load in current filt sng and .not.mat
    crntfile = fscanf(fidB,'%s',1)
    if isempty(crntfile)
        break
    end
    filtfile = [char(crntfile) '.filt'];
    [filtsong,fs]=read_filt(filtfile);
    notmatfile = [char(crntfile) '.not.mat'];
    load(char(notmatfile));
    %make sure to get all sound
    ons = fix((onsets-5).*fs/1000);
    offs = fix((offsets+5).*fs/1000);
    % ensure ids are with in filtsong
    ons(ons<=0) = 1; offs(offs>=length(filtsong)) = length(filtsong);

    % calculate feature vectors for appropriate labels
    featvect = [];
    for i = 1:size(ons,1)
        if ~isempty(strfind(char(labelIndex)',char(labels(i))))
            featvect = [featvect; feature_vect(filtsong(ons(i):offs(i)),fs)];
        else
            featvect = [featvect; zeros(1,length(template_mean(1,:)))];
        end
    end

    %now find distances
    for j = 1:length(labels)
        if sum(featvect(j,:)) ~= 0
            k = strfind(char(labelIndex'),labels(j));
            kmean=template_mean(k,:);
            kdiff=featvect(j,:)-kmean;
            kdiff(:,6)=kdiff(:,6)+100*(kdiff(:,6)<-50)-100*(kdiff(:,6)>50);
            kdiff(:,7)=kdiff(:,7)+100*(kdiff(:,7)<-50)-100*(kdiff(:,7)>50);
            kinvcov=inv(template_cov{k});
            dist_crnt = sqrt(sum((kdiff*kinvcov).*kdiff,2));
            if min==1
                if dist_crnt < d(k).top10(end) % store the minimum 10 syls
                    d(k).top10(end) = dist_crnt;
                    mahalmin{k}.syl{h}.snd = filtsong(ons(j):offs(j));
                    mahalmin{k}.dist(end) = dist_crnt;
                    [d(k).top10,idx] = sort(d(k).top10);
                    mahalmin{k}.label = labelIndex(k);
                    mahalmin{k}.FS = fs;
                    for h = 1:10
                        l = idx(h);
                        crntsnd = mahalmin{k}.syl{l}.snd;
                        mahalmin{k}.syl{h}.snd = crntsnd;
                        mahalmin{k}.dist(h) = mahalmin{k}.dist(l);
                    end
                end
            else
                if dist_crnt < d(k).top10(end) % store the minimum 10 syls
                    d(k).top10(end) = dist_crnt;
                    mahalmin{k}.syl{h}.snd = filtsong(ons(j):offs(j));
                    mahalmin{k}.dist(end) = dist_crnt;
                    [d(k).top10,idx] = sort(d(k).top10);
                    mahalmin{k}.label = labelIndex(k);
                    mahalmin{k}.FS = fs;
                    for h = 1:10
                        l = idx(h);
                        crntsnd = mahalmin{k}.syl{l}.snd;
                        mahalmin{k}.syl{h}.snd = crntsnd;
                        mahalmin{k}.dist(h) = mahalmin{k}.dist(l);
                    end
                end
            end
        end
    end
end

cnt2 = 0;
for r = 1:length(labelIndex)
    if ~isempty(mahalmin{r}.syl{1}.snd)
        cnt2 = cnt2+1;
        mahalmin_fin{cnt2} = mahalmin{r};
    end
end

if min==1
    name2save = [char(name) '.mhlsyl.min_strct.mat'];
else
    name2save = [char(name) '.mhlsyl.max_strct.mat'];
end
save(char(name2save), 'mahalmin_fin');

make_mahal_stim(mahalmin_fin,name2save);

