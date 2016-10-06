function f = psygumbelcdf(x,mu,sigma)
%PSYGUMBELCDF Cumulative Gumbel distribution psychometric curve.

% Convert parametrization from median and standard deviation
beta = sqrt(6)/pi*sigma;
mu = bsxfun(@plus, mu, beta*(-0.366512920581664)); % log(log(2))

z = bsxfun(@rdivide,bsxfun(@minus,x,mu),beta);
f = exp(-exp(-z));

end