function [xnext,psy,output] = psybayes_joint(psy,method,vars,xi,yi,ci)
%PSYBAYES_JOINT Joint Bayesian adaptive estimation of psychometric functions.
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
if nargin < 5; ci = []; end

persistent firstcall;

xnext = [];

if isempty(firstcall)
    firstcall = 0;        
    % Add all subdirectories to MATLAB path
    [path,~,~] = fileparts(mfilename('fullpath'));
    addpath(genpath(path));
end

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


%% Initialization of PSY structures

% Empty struct, a single psychometric function
if isempty(psy); psy = {[]}; end

% PSY can be NCND, number of experimental conditions
if isnumeric(psy) && isscalar(psy)
    psy = cell(1,psy);    
end

Ncnd = numel(psy);          % Number of experimental conditions
sharedlambda = Ncnd > 1;    % If fitting multiple conditions, assume lapse rate is shared

% Number of psychometric curves
Nfuns = zeros(1,Ncnd);

for c = 1:Ncnd
    if isempty(psy{c}) || ~isfield(psy{c},'post')

        % Call initialization function
        psyinfo = psy{c};
        [psy{c},Nfuns(c)] = psyinit(psyinfo,Ncnd);

        % Enforce symmetry of test stimuli? (symmetric wrt left/right of the 
        % mean of the psychometric curve)
        if isfield(psyinfo,'forcesymmetry') && ~isempty(psyinfo.forcesymmetry)
            psy{c}.forcesymmetry = psyinfo.forcesymmetry;
        else
            psy{c}.forcesymmetry = 0;  
        end    
    else
        % Reset psychometric function
        [psy{c},Nfuns(c)] = psyfunset(psy{c});
    end
end

if ~all(Nfuns == Nfuns(1))
    error('All conditions should have the same number of psychometric curves.');
end
Nfuns = Nfuns(1);

if Ncnd > 1 && Nfuns > 1
    error('For the moment joint psychometric curve fitting only supports a single psychometric curve model.');
end

% Initialize psychometric functions in each condition
for c = 1:Ncnd        
    % Convert psychometric function to function handle
    if ~iscell(psy{c}.psychofun)
        psychofun{c}{1} = str2func(psy{c}.psychofun);
    else
        for k = 1:Nfuns
            psychofun{c}{k} = str2func(psy{c}.psychofun{k});
        end
    end
    
    % Precompute psychometric function
    if isempty(psy{c}.f)
        for k = 1:Nfuns
            psy{c}.f{k} = psychofun{c}{k}(psy{c}.x,psy{c}.mu,psy{c}.sigma,psy{c}.lambda,psy{c}.gamma);
            % Check if last stimulus is easy stimulus (by default Inf)
            if psy{c}.x(end) == Inf
                if isempty(psy{c}.gamma)
                    temp(1,1,:) = 1-psy{c}.lambda/2;            
                else
                    temp(1,1,:) = 1-psy{c}.lambda*(1-psy{c}.gamma);
                end
                psy{c}.f{k}(:,:,:,end) = repmat(temp,[numel(psy{c}.mu),numel(psy{c}.logsigma),1]);
            end
        end
    end    
end

% Update log posterior given the new data points XI, YI
if ~isempty(xi) && ~isempty(yi)
    
    if isempty(ci) && Ncnd == 1
        ci = 1;
    elseif isempty(ci)
        error('Current condition index CI not specified.');
    elseif ~isscalar(ci)
        error('Current condition index needs to be a scalar.');
    end
    
    for k = 1:Nfuns
        for i = 1:numel(xi)            
            cii = ci;
            
            % Maximum precision stimulus
            if isinf(xi(i))
                if isempty(psy{cii}.gamma)
                    if yi(i) == 1
                        like = 1-psy{cii}.lambda/2;
                    elseif yi(i) == 0
                        like = psy{cii}.lambda/2;
                    end                    
                else
                    if yi(i) == 1
                        like = 1-psy{cii}.lambda*(1-psy{cii}.gamma);
                    elseif yi(i) == 0
                        like = psy{cii}.lambda*(1-psy{cii}.gamma);
                    end
                end
            else
                if yi(i) == 1
                    like = psychofun{cii}{k}(xi(i),psy{cii}.mu,psy{cii}.sigma,psy{cii}.lambda,psy{cii}.gamma);
                elseif yi(i) == 0
                    like = 1 - psychofun{cii}{k}(xi(i),psy{cii}.mu,psy{cii}.sigma,psy{cii}.lambda,psy{cii}.gamma);
                end
            end
            
            % Save unnormalized log posterior
            psy{cii}.logupost{k} = bsxfun(@plus, psy{cii}.logupost{k}, log(like));
            
            % Compute posterior over lambda for this condition
            temp = psy{cii}.logupost{k};
            temp = exp(temp - max(temp(:)));
            psy{cii}.postlambda{k} = sum(sum(temp,1),2) / sum(temp(:));
            
            % Compute joint posterior over lambda
            if sharedlambda
                postlambda_joint{k} = ones(size(psy{cii}.postlambda{k}));
                for c = 1:Ncnd
                    postlambda_joint{k} = postlambda_joint{k} .* psy{c}.postlambda{k};
                end
                postlambda_joint{k} = postlambda_joint{k} / sum(postlambda_joint{k});
            end

            % Compute normalized posterior
            temp = exp(psy{cii}.logupost{k} - max(psy{cii}.logupost{k}(:)));
            if sharedlambda
                temp = bsxfun(@times, temp, postlambda_joint{k} ./ psy{cii}.postlambda{k});
            end
            psy{cii}.post{k} = temp./sum(temp(:));
        end
    end
    
    % Update data
    for i = 1:numel(xi)
        cii = ci;
        psy{cii}.ntrial = psy{cii}.ntrial + 1;
        psy{cii}.data = [psy{cii}.data; xi(i) yi(i)];
    end
        
    % Update refractory times list for each presented stimulus
    for i = 1:numel(xi)
        cii = ci;
        psy{cii}.reflist = max(psy{cii}.reflist - 1, 0);
        if psy{cii}.reftime > 0 && isfinite(xi(i))
            idx = (psy{cii}.x <= xi(i) + psy{cii}.refradius) & (psy{cii}.x >= xi(i) - psy{cii}.refradius);
            wtrials(1,1,1,:) = geornd(1/(1+psy{cii}.reftime)*ones(1,sum(idx)));
            psy{cii}.reflist(idx) = max(wtrials, psy{cii}.reflist(idx));
        end
    end    
end

% Compute posterior over psychometric functions
for c = 1:Ncnd
    if Nfuns > 1
        logp = zeros(1,Nfuns);
        for k = 1:Nfuns
            logp(k) = logsumexp(psy{c}.logupost{k}(:));
        end        
        % This is not correct for shared lambda
        psy{c}.psychopost = exp(logp - max(logp));
        psy{c}.psychopost = psy{c}.psychopost ./ sum(psy{c}.psychopost);
    else
        psy{c}.psychopost = 1;
    end
end

% Only one argument assumes that this is the final call
if nargin < 2
    for c = 1:Ncnd
        psy{c}.f = [];     % Empty some memory
        psy{c}.reflist = zeros(size(psy{c}.x));   % Reset refractory times list
    end
    return;
end

% Compute mean of the posterior of mu for the current condition
postmu = zeros(numel(psy{ci}.mu),Nfuns);
for k = 1:Nfuns
    postmu(:,k) = sum(sum(psy{ci}.post{k},2),3);
end
emu = sum(sum(bsxfun(@times, bsxfun(@times, psy{ci}.psychopost, postmu), psy{ci}.mu),2),1);

% Randomly remove half of the x
if psy{ci}.forcesymmetry
    if rand() < 0.5; xindex = psy{ci}.x < emu; else xindex = psy{ci}.x >= emu; end
else
    xindex = true(size(psy{ci}.x));
end

% Consider only available stimuli
xindex = xindex & (psy{ci}.reflist == 0);

% No stimuli are available, free some stimuli and reset refractory list
if all(xindex == 0)
    xindex(psy{ci}.reflist == min(psy{ci}.reflist)) = 1;
    psy{ci}.reflist = zeros(size(psy{ci}.x));
end

% Compute sampling point X that minimizes expected chosen criterion
if nargin > 0
        
    Nx = numel(psy{ci}.x);
        
    r1 = zeros(1,1,1,Nx,Nfuns);
    post1 = zeros([size(psy{ci}.post{1}),Nx,Nfuns]);
    post0 = zeros([size(psy{ci}.post{1}),Nx,Nfuns]);
    if Nfuns > 1
        u1 = zeros(1,1,1,Nx,Nfuns);
        u0 = zeros(1,1,1,Nx,Nfuns);
    end
    
    for k = 1:Nfuns
        if Nfuns > 1
            % Compute posteriors and unnormalized model evidence at next step for R=1 and R=0
%            [post1(:,:,:,:,k),post0(:,:,:,:,k),r1(1,1,1,:,k),u1(1,1,1,:,k),u0(1,1,1,:,k)] = nextposterior(psy.f{k}(:,:,:,xindex),psy.post{k},psy.logupost{k});
            % This is not completely correct for shared lapse rate
            [post1(:,:,:,:,k),post0(:,:,:,:,k),r1(1,1,1,:,k),u1(1,1,1,:,k),u0(1,1,1,:,k)] = nextposterior(psy{ci}.f{k}(:,:,:,:),psy{ci}.post{k},psy{ci}.logupost{k});
        else
            % Compute posteriors at next step for R=1 and R=0
            %[post1(:,:,:,:,k),post0(:,:,:,:,k),r1(1,1,1,:,k)] = nextposterior(psy.f{k}(:,:,:,xindex),psy.post{k});
            [post1(:,:,:,:,k),post0(:,:,:,:,k),r1(1,1,1,:,k)] = nextposterior(psy{ci}.f{k}(:,:,:,:),psy{ci}.post{k});
        end
    end
    
    % Marginalize over unrequested variables
    index = find(~vars);
    for iTheta = index
        post1 = sum(post1,iTheta);
        post0 = sum(post0,iTheta);
    end
        
    if Nfuns > 1
        u0 = exp(bsxfun(@minus, u0, max(u0,[],5)));
        u0 = bsxfun(@rdivide,u0,sum(u0,5));
        u1 = exp(bsxfun(@minus, u1, max(u1,[],5)));
        u1 = bsxfun(@rdivide,u1,sum(u1,5));
        post1 = bsxfun(@times, u1, post1);
        post0 = bsxfun(@times, u0, post0);
        w(1,1,1,1,:) = psy{ci}.psychopost;
        r1 = sum(bsxfun(@times, w, r1), 5);
    end
    
    switch lower(method)
        case {'var','variance'}            
            post0 = squeeze(sum(post0, 5));
            post1 = squeeze(sum(post1, 5));            
            index = find(vars,1);
            switch index
                case 1; qq = psy{ci}.mu(:);
                case 2; qq = psy{ci}.logsigma(:);
                case 3; qq = psy{ci}.lambda(:);
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
            
        case {'proj','projection'}
            Nsteps = 4; 
            max(r1)
            anchors = linspace(0.5, max(r1), Nsteps+1);
            anchors = anchors(2:end);
            target = zeros(Nx,1);
            for jj = 1:numel(anchors)
                for k = 1:Nfuns
                    [~,idx(:,:,:,1,k)] = min(abs(psy{ci}.f{k}-anchors(jj)),[],4);
                end
                mean1 = sum(sum(sum(sum(bsxfun(@times,post1,idx),1),2),3),5);
                mean0 = sum(sum(sum(sum(bsxfun(@times,post0,idx),1),2),3),5);
                var1 = sum(sum(sum(sum(bsxfun(@times,post1,idx.^2),1),2),3),5) - mean1.^2;
                var0 = sum(sum(sum(sum(bsxfun(@times,post0,idx.^2),1),2),3),5) - mean0.^2;
                target = target + r1(:).*var1(:) + (1-r1(:)).*var0(:);
            end
            
        case {'model'}
            H1 = -u1.*log(u1);
            H0 = -u0.*log(u0);            
            H1(~isfinite(H1)) = 0;
            H0(~isfinite(H0)) = 0;
            H1 = sum(H1,5);
            H0 = sum(H0,5);
            target = r1(:).*H1(:) + (1-r1(:)).*H0(:);
            
            
        otherwise
            error('Unknown method. Allowed methods are ''var'' and ''ent'' for, respectively, predicted variance and predicted entropy minimization.');
    end

    % Store target for plotting
    psy{ci}.target = target(:)';
    
    % Location X that minimizes target metric
    [~,index] = min(target(xindex));
    xred = psy{ci}.x(xindex);
    xnext = xred(index);
    
    psy{ci}.xnext = xnext;
end

% Compute parameter estimates
if nargout > 2
    w = psy{ci}.psychopost;
    
    % Compute mean and variance of the estimate of MU
    postmu = marginalpost(psy{ci}.post,w,[2,3]);
    postmu = postmu./sum(postmu,1);
    emu = sum(postmu.*psy{ci}.mu,1);
    estd = sqrt(sum(postmu.*psy{ci}.mu.^2,1) - emu.^2);
    output.mu.mean = emu;
    output.mu.std = estd;
    
    % Compute mean and variance of the estimate of LOGSIGMA and SIGMA
    postlogsigma = marginalpost(psy{ci}.post,w,[1,3]);
    postlogsigma = postlogsigma./sum(postlogsigma,2);    
    emu = sum(postlogsigma.*psy{ci}.logsigma,2);
    estd = sqrt(sum(postlogsigma.*psy{ci}.logsigma.^2,2) - emu.^2);
    output.logsigma.mean = emu;
    output.logsigma.std = estd;
    
    postsigma = postlogsigma./psy{ci}.sigma;
    postsigma = postsigma./sum(postsigma,2);    
    emu = sum(postsigma.*psy{ci}.sigma,2);
    estd = sqrt(sum(postsigma.*psy{ci}.sigma.^2,2) - emu.^2);
    output.sigma.mean = emu;
    output.sigma.std = estd;
    
    % Compute mean and variance of the estimate of LAMBDA
    if sharedlambda && Nfuns == 1
        postlambda = psy{ci}.postlambda{1};        
    else
        % This is possibly not correct
        postlambda = marginalpost(psy{ci}.post,w,[1,2]);
    end
    postlambda = postlambda./sum(postlambda,3);
    emu = sum(postlambda.*psy{ci}.lambda,3);
    estd = sqrt(sum(postlambda.*psy{ci}.lambda.^2,3) - emu.^2);
    output.lambda.mean = emu;
    output.lambda.std = estd;
end

end


%--------------------------------------------------------------------------
function [post1,post0,r1,u1,u0] = nextposterior(f,post,logupost)
%NEXTPOSTERIOR Compute posteriors on next trial depending on possible outcomes

    mf = 1-f;
    post1 = bsxfun(@times, post, f);
    r1 = sum(sum(sum(post1,1),2),3);
    post0 = bsxfun(@times, post, mf);
    post1 = bsxfun(@rdivide, post1, sum(sum(sum(post1,1),2),3));
    post0 = bsxfun(@rdivide, post0, sum(sum(sum(post0,1),2),3));    
    
    if nargin > 2 && nargout > 3
        logupost1 = bsxfun(@plus, logupost, log(f));
        logupost0 = bsxfun(@plus, logupost, log(mf));
        z0 = max(logupost0(:));
        u0 = log(sum(sum(sum(exp(logupost0 - z0),1),2),3));
        z1 = max(logupost1(:));
        u1 = log(sum(sum(sum(exp(logupost1 - z1),1),2),3));        
    end
end

%--------------------------------------------------------------------------
function r = geornd(p)
%GEORND Random arrays from the geometric distribution.

p(p <= 0 | p > 1) = NaN;    % Return NaN for illegal parameter values
r = ceil(abs(log(rand(size(p))) ./ log(1 - p)) - 1); % == geoinv(u,p)
r(r < 0) = 0;   % Force a zero when p==1, instead of -1

end