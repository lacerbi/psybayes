function [xnext,psy,output] = psybayes(psy,method,vars,xi,yi)
%PSYBAYES Bayesian adaptive estimation of psychometric function.
%
%  PSYBAYES implements Kontsevich and Tyler's (1999) Bayesian adaptive 
%  method PSI for estimation of parameters of the psychometric function via
%  maximization of information gain (including lapse; see Prins 2012). 
%  PSYBAYES also supports the marginal-PSI method by Prins (2013).
%
%  See PSYTEST for documentation and a working usage example.
%
%  References: 
%  Kontsevich, L. L., & Tyler, C. W. (1999). "Bayesian adaptive estimation 
%  of psychometric slope and threshold". Vision Research, 39(16), 2729-2737.
%
%  Prins, N. (2012). "The adaptive psi method and the lapse rate". Journal 
%  of Vision, 12(9), 322-322. (link)
%
%  Prins, N. (2013). "The psi-marginal adaptive method: How to give nuisance 
%  parameters the attention they deserve (no more, no less)". Journal of 
%  Vision, 13(7), 3-3.
%
%  See also PSYBAYES_PLOT, PSYTEST.

% Copyright (C) 2016 Luigi Acerbi
%
% This software is distributed under the GNU General Public License 
% (version 3 or later); please refer to the file LICENSE.txt, included with 
% the software, for details.

%   Author:     Luigi Acerbi
%   Email:      luigi.acerbi@gmail.com
%   Version:    24/Jul/2016


if nargin < 1; psy = []; end
if nargin < 2; method = []; end
if nargin < 3; vars = []; end
if nargin < 4; xi = []; yi = []; end

xnext = [];

% Default method is expected entropy minimization
if isempty(method); method = 'ent'; end

if isempty(vars)
    switch lower(method)
        case 'ent'; vars = [1 1 1];
        case 'var'; vars = [1 0 0];
        otherwise
            error('Unknown optimization method.');
    end
end
if numel(vars) ~= 3; error('VARS need to be a 3-element array for MU, SIGMA and LAMBDA.'); end

%% First call, initialize everything
if isempty(psy) || ~isfield(psy,'post')
    
    % Call initialization function
    psyinfo = psy;
    psy = psyinit(psyinfo);
    
    % Enforce symmetry of test stimuli? (symmetric wrt left/right of the 
    % mean of the psychometric curve)
    if isfield(psyinfo,'forcesymmetry') && ~isempty(psyinfo.forcesymmetry)
        psy.forcesymmetry = psyinfo.forcesymmetry;
    else
        psy.forcesymmetry = 0;  
    end    
 
end

% Choose correct psychometric function (YES/NO or PCORRECT)
if ~isempty(psy.gamma)
    psychofun = @(x_,mu_,sigma_,lambda_) psyfun_pcorrect(x_,mu_,sigma_,lambda_,psy.gamma);
else    
    psychofun = @(x_,mu_,sigma_,lambda_) psyfun_yesno(x_,mu_,sigma_,lambda_);
end

% Precompute psychometric function
if isempty(psy.f) || isempty(psy.mf)
    psy.f = psychofun(psy.x,psy.mu,psy.sigma,psy.lambda);
    psy.mf = 1 - psy.f;
end

% Update log posterior given the new data points XI, YI
if ~isempty(xi) && ~isempty(yi)    
    for i = 1:numel(xi)
        if yi(i) == 1
            like = psychofun(xi(i),psy.mu,psy.sigma,psy.lambda);
        elseif yi(i) == 0
            like = 1 - psychofun(xi(i),psy.mu,psy.sigma,psy.lambda); 
        end
        psy.post = psy.post.*like;
    end
    psy.post = psy.post./sum(psy.post(:));
    
    psy.ntrial = psy.ntrial + numel(xi);
    psy.data = [psy.data; xi(:) yi(:)];
end

% Compute mean of the posterior of mu
postmu = sum(sum(psy.post,2),3);
emu = sum(postmu.*psy.mu);

% Randomly remove half of the x
if psy.forcesymmetry
    if rand() < 0.5; xindex = psy.x < emu; else xindex = psy.x >= emu; end
else
    xindex = true(size(psy.x));
end

% Compute sampling point X that minimizes expected chosen criterion
if nargin > 0
        
    xred = psy.x(xindex);
    
    % Compute posteriors at next step for R=1 and R=0
    [post1,post0,r1] = nextposterior(psy.f(:,:,:,xindex),psy.post);

    % Marginalize over unrequested variables
    index = find(~vars);
    for iTheta = index
        post1 = sum(post1,iTheta);
        post0 = sum(post0,iTheta);
    end
    
    switch lower(method)
        case {'var','variance'}
            post1 = squeeze(post1);
            post0 = squeeze(post0);
            index = find(vars,1);
            switch index
                case 1; qq = psy.mu(:);
                case 2; qq = psy.logsigma(:);
                case 3; qq = psy.lambda(:);
            end
            mean1 = sum(bsxfun(@times,post1,qq),1);
            mean0 = sum(bsxfun(@times,post0,qq),1);
            var1 = sum(bsxfun(@times,post1,qq.^2),1) - mean1.^2;
            var0 = sum(bsxfun(@times,post0,qq.^2),1) - mean0.^2;
            target = r1(:).*var1(:) + (1-r1(:)).*var0(:);
                        
        case {'ent','entropy'}
            temp1 = -post1.*log(post1);
            temp0 = -post0.*log(post0);            
            temp1(~isfinite(temp1)) = 0;
            temp0(~isfinite(temp0)) = 0;
            H1 = temp1;     H0 = temp0;
            for iTheta = find(vars)
                H1 = sum(H1,iTheta);
                H0 = sum(H0,iTheta);
            end
            target = r1(:).*H1(:) + (1-r1(:)).*H0(:);
            
        otherwise
            error('Unknown method. Allowed methods are ''var'' and ''ent'' for, respectively, predicted variance and predicted entropy minimization.');
    end

    % Location X that minimizes target metric
    [~,index] = min(target);
    xnext = xred(index);
    
    psy.xnext = xnext;
end

% Compute parameter estimates
if nargout > 2
    
    % Compute mean and variance of the estimate of MU
    postmu = sum(sum(psy.post,2),3);
    postmu = postmu./sum(postmu,1);
    emu = sum(postmu.*psy.mu,1);
    estd = sqrt(sum(postmu.*psy.mu.^2,1) - emu.^2);
    output.mu.mean = emu;
    output.mu.std = estd;
    
    % Compute mean and variance of the estimate of LOGSIGMA and SIGMA
    postlogsigma = sum(sum(psy.post,1),3);    
    postlogsigma = postlogsigma./sum(postlogsigma,2);    
    emu = sum(postlogsigma.*psy.logsigma,2);
    estd = sqrt(sum(postlogsigma.*psy.logsigma.^2,2) - emu.^2);
    output.logsigma.mean = emu;
    output.logsigma.std = estd;
    
    postsigma = postlogsigma./psy.sigma;
    postsigma = postsigma./sum(postsigma,2);    
    emu = sum(postsigma.*psy.sigma,2);
    estd = sqrt(sum(postsigma.*psy.sigma.^2,2) - emu.^2);
    output.sigma.mean = emu;
    output.sigma.std = estd;
    
    % Compute mean and variance of the estimate of LAMBDA
    postlambda = sum(sum(psy.post,1),2);
    postlambda = postlambda./sum(postlambda,3);    
    emu = sum(postlambda.*psy.lambda,3);
    estd = sqrt(sum(postlambda.*psy.lambda.^2,3) - emu.^2);
    output.lambda.mean = emu;
    output.lambda.std = estd;
end

% Only one argument assumes that this is the final call
if nargin < 2
    % Empty some memory
    psy.f = []; psy.mf = [];
end

end


%--------------------------------------------------------------------------
function [post1,post0,r1] = nextposterior(f,post)
%NEXTPOSTERIOR Compute posteriors on next trial depending on possible outcomes

    mf = 1-f;
    post1 = bsxfun(@times, post, f);
    r1 = sum(sum(sum(post1,1),2),3);
    post0 = bsxfun(@times, post, mf);
    post1 = bsxfun(@rdivide, post1, sum(sum(sum(post1,1),2),3));
    post0 = bsxfun(@rdivide, post0, sum(sum(sum(post0,1),2),3));    
end