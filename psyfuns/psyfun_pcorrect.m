function f = psyfun_pcorrect(x,mu,sigma,lambda,gamma,F)
%PSYFUN_PCORRECT Psychometric function for percent correct (with guessing)

f = bsxfun(@plus, gamma, ...
    bsxfun(@times,1-gamma-lambda+gamma*lambda, F(x,mu,sigma)));

end