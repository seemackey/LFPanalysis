function [trigtype,trigtimes] = EphysExtractTrigs(trig,trigch,newadrate)

    % noise - 1, led (open eyes) - 3, sacc on 4, sacc off 5
    trig0=trig.anatrig{trigch};


    %triggers if there are multiple, trigch is the one we're selecting


    if trigch==1  % aud stim is stored as cell type SOMETIMES so we need a lot of ifs

        if iscell(trig.ttype)
            present_trigidx=trig.ttype{1,trigch}(:,:)>0; % places where there was a trigger
            trigtype=trig.ttype{1,trigch}(present_trigidx); % aud stim are in cells usually
        else
            present_trigidx=trig.ttype(:,:)>0;
            trigtype=trig.ttype(present_trigidx); 
        end

    elseif trigch==2
        present_trigidx=trig.ttype(:,:)>0;
        trigtype=trig.ttype(present_trigidx); % vis stim are in double ... should be fixed

    elseif trigch==3 % led stim when eyes are open, gotta do something special here

        trig0_all=trig.anatrig{2}; %ALL (eyes open & close) trigs
        for i=1:length(trig0)
           trig_mutual(i)= find(trig0_all==trig0(i)); %gives you the ALL trig number that is in trig0
        end
        trigsorttmp=trig.triglength{1,2};
        trigtype(:,1)=trigsorttmp(trig_mutual); % std/dev when eyes are open % vis stim are in double ... should be fixed


    elseif trigch==4 % sacc on
        trigtype=ones(size(trig.anatrig{1,trigch}(:)));% just sticking ones here because we don't have diff sacc types

    elseif trigch==5 % sacc off
        trigtype=ones(size(trig.anatrig{1,trigch}(:))); % just sticking ones here because we don't have diff sacc types
    end

    trigtimes = [];
    for trigredxct=1:length(trig0)
     trigtimes(trigredxct)    = round(trig0(trigredxct)./(trig.adrate/newadrate));
    end

end