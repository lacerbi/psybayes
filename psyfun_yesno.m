function f = psyfun_yesno(x,mu,sigma,lambda)
%PSYFUN_YESNO Psychometric function for YES/NO tasks

f = bsxfun(@plus, lambda/2, ...
    bsxfun(@times,1-lambda,0.5*(1+erf(bsxfun(@rdivide,bsxfun(@minus,x,mu),sqrt(2)*sigma)))));

end