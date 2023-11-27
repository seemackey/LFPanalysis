function [norm] = jigglenorm(data)
%% give data in single column
%% assumes you're normalizing row 1:4 to row 1

    for normct = 1:4:length(data)

        norm(normct:normct+3,1) = data(normct:normct+3,1)/data(normct,1);
        % if norm(normct+1:normct+3) is too big or too small we delete

    end


end