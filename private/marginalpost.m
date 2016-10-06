function y = marginalpost(post,w,idx)
%MARGINALPOST Compute marginal posterior

    Nfuns = numel(post);
    for k = 1:Nfuns
        for j = idx
            post{k} = sum(post{k},j);
        end
    end    
    y = zeros(size(post{1}));
    for k = 1:Nfuns; y = y + w(k)*post{k}; end
end