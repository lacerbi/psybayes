function f = psyfun_pcorrect(x,mu,sigma,lambda,gamma)
%PSYFUN_PCORRECT Psychometric function for percent correct (with guessing)

f = bsxfun(@plus, gamma, ...
    bsxfun(@times,1-gamma-lambda,0.5*(1+erf(bsxfun(@rdivide,bsxfun(@minus,x,mu),sqrt(2)*sigma)))));

end