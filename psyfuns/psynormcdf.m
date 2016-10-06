function f = psynormcdf(x,mu,sigma)
%PSYNORMCDF Cumulative normal distribution psychometric curve.

z = bsxfun(@rdivide,bsxfun(@minus,x,mu),sqrt(2)*sigma);
f = 0.5 + 0.5*erf(z);

end