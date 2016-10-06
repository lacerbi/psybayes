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
%   Version:    05/Oct/2016

if nargin < 1; psy = []; end
if nargin < 2; method = []; end
if nargin < 3; vars = []; end
if nargin < 4; xi = []; yi = []; end

xnext = [];

% Default method is expected entropy minimization
if isempty(method); method = 'ent'; end

% Marginal-PSI method, select parameters of interest
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
    [psy,Nfuns] = psyinit(psyinfo);
    
    % Enforce symmetry of test stimuli? (symmetric wrt left/right of the 
    % mean of the psychometric curve)
    if isfield(psyinfo,'forcesymmetry') && ~isempty(psyinfo.forcesymmetry)
        psy.forcesymmetry = psyinfo.forcesymmetry;
    else
        psy.forcesymmetry = 0;  
    end    
else
    % Reset psychometric function
    [psy,Nfuns] = psyfunset(psy);
end

% Select psychometric function
if ~iscell(psy.psychofun)
    psychofun{1} = str2func(psy.psychofun);
else
    for k = 1:Nfuns
        psychofun{k} = str2func(psy.psychofun{k});
    end
end
    
% Precompute psychometric function
if isempty(psy.f) || isempty(psy.mf)
    for k = 1:Nfuns
        psy.f{k} = psychofun{k}(psy.x,psy.mu,psy.sigma,psy.lambda,psy.gamma);
        psy.mf{k} = 1 - psy.f{k};
    end
end

% Update log posterior given the new data points XI, YI
if ~isempty(xi) && ~isempty(yi)
    for k = 1:Nfuns
        for i = 1:numel(xi)
            if isinf(xi(i))
                if isempty(psy.gamma)
                    if yi(i) == 1
                        like = 1-psy.lambda/2;
                    elseif yi(i) == 0
                        like = psy.lambda/2;
                    end                    
                else
                    if yi(i) == 1
                        like = 1-psy.lambda*(1-psy.gamma);
                    elseif yi(i) == 0
                        like = psy.lambda*(1-psy.gamma);
                    end
                end
            else
                if yi(i) == 1
                    like = psychofun{k}(xi(i),psy.mu,psy.sigma,psy.lambda,psy.gamma);
                elseif yi(i) == 0
                    like = 1 - psychofun{k}(xi(i),psy.mu,psy.sigma,psy.lambda,psy.gamma);
                end
            end
            
            % Save unnormalized log posterior
            psy.logupost{k} = bsxfun(@plus, psy.logupost{k}, log(like));

            % Compute normalized posterior
            psy.post{k} = exp(psy.logupost{k} - max(psy.logupost{k}(:)));            
            psy.post{k} = psy.post{k}./sum(psy.post{k}(:));
        end
    end
    
    psy.ntrial = psy.ntrial + numel(xi);
    psy.data = [psy.data; xi(:) yi(:)];
end

% Compute posterior over psychometric functions
if Nfuns > 1
    logp = log(psy.psychoprior);
    for k = 1:Nfuns
        logp(k) = logp(k) + logsumexp(psy.logupost{k}(:));
    end
    psy.psychopost = exp(logp - max(logp));
    psy.psychopost = psy.psychopost ./ sum(psy.psychopost);
    psy.psychopost
else
    psy.psychopost = 1;
end

% Compute mean of the posterior of mu
postmu = zeros(numel(psy.mu),Nfuns);
for k = 1:Nfuns
    postmu(:,k) = sum(sum(psy.post{k},2),3);
end
emu = sum(sum(bsxfun(@times, bsxfun(@times, psy.psychopost, postmu), psy.mu),2),1);

% Randomly remove half of the x
if psy.forcesymmetry
    if rand() < 0.5; xindex = psy.x < emu; else xindex = psy.x >= emu; end
else
    xindex = true(size(psy.x));
end

% Compute sampling point X that minimizes expected chosen criterion
if nargin > 0
    
    Nx = numel(psy.x);
        
    xred = psy.x(xindex);
    r1 = zeros(1,1,1,Nx,Nfuns);
    post1 = zeros([size(psy.post{1}),Nx,Nfuns]);
    post0 = zeros([size(psy.post{1}),Nx,Nfuns]);
    
    for k = 1:Nfuns
        % Compute posteriors at next step for R=1 and R=0
        [post1(:,:,:,:,k),post0(:,:,:,:,k),r1(1,1,1,:,k)] = nextposterior(psy.f{k}(:,:,:,xindex),psy.post{k});
    end
    
    % Marginalize over unrequested variables
    index = find(~vars);
    for iTheta = index
        post1 = sum(post1,iTheta);
        post0 = sum(post0,iTheta);
    end    
        
    if Nfuns > 1
        w(1,1,1,1,:) = psy.psychopost;
        post1 = bsxfun(@times, w, post1);
        post0 = bsxfun(@times, w, post0);
        r1 = sum(bsxfun(@times, w, r1), 5);
    end
    
    switch lower(method)
        case {'var','variance'}            
            post0 = sum(bsxfun(@times, post0, w), 5);
            post1 = sum(bsxfun(@times, post1, w), 5);            
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
            if Nfuns > 1
                H1 = sum(H1,5);
                H0 = sum(H0,5);
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
    w = psy.psychopost;
    
    % Compute mean and variance of the estimate of MU
    postmu = marginalpost(psy.post,w,[2,3]);
    postmu = postmu./sum(postmu,1);
    emu = sum(postmu.*psy.mu,1);
    estd = sqrt(sum(postmu.*psy.mu.^2,1) - emu.^2);
    output.mu.mean = emu;
    output.mu.std = estd;
    
    % Compute mean and variance of the estimate of LOGSIGMA and SIGMA
    postlogsigma = marginalpost(psy.post,w,[1,3]);
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
    postlambda = marginalpost(psy.post,w,[1,2]);
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