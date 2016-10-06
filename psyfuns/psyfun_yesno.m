function f = psyfun_yesno(x,mu,sigma,lambda,gamma,F)
%PSYFUN_YESNO Psychometric function for YES/NO tasks

f = bsxfun(@plus, lambda/2, ...
    bsxfun(@times,1-lambda, F(x,mu,sigma)));

end