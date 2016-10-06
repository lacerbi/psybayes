function f = psylogicdf(x,mu,sigma)
%PSYLOGICDF Cumulative logistic distribution psychometric curve.

s = sqrt(3)/pi*sigma;
z = bsxfun(@rdivide,bsxfun(@minus,x,mu),s);
f = 1./(1+exp(-z));

end