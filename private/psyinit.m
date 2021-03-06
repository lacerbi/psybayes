function [psy,Nfuns] = psyinit(psyinfo,Ncnd)
%PSYINIT Initialize PSY struct.

% Total number of conditions (one by default)
if nargin < 2 || isempty(Ncnd); Ncnd = 1; end

psy = [];

psy.ntrial = 0;     % Trial number
psy.data = [];      % Record of data

if ~isfield(psyinfo,'psychofun'); psyinfo.psychofun = []; end
if iscell(psyinfo.psychofun)
    Nfuns = numel(psyinfo.psychofun);
    cellflag = 1;
else
    Nfuns = 1;
    cellflag = 0;
end

% Default grid sizes
K = 1/(Nfuns^0.25);
nx = round(65*K);
nmu = round(51*K);
nsigma = round(25*K);
nlambda = round(25*K);

psy.mu = [];
psy.logsigma = [];
psy.lambda = [];
psy.x = [];

% Grid over parameters of psychometric function
if isfield(psyinfo,'range')

    % Get grid sizes (third element in initialization range field)
    if isfield(psyinfo.range,'mu') && numel(psyinfo.range.mu > 2)
        nmu = psyinfo.range.mu(3);
    end
    if isfield(psyinfo.range,'sigma') && numel(psyinfo.range.sigma > 2)
        nsigma = psyinfo.range.sigma(3);
    elseif isfield(psyinfo.range,'logsigma') && numel(psyinfo.range.logsigma > 2)
        nsigma = psyinfo.range.logsigma(3);            
    end
    if isfield(psyinfo.range,'lambda') && numel(psyinfo.range.lambda > 2)
        nlambda = psyinfo.range.lambda(3);
    end
    if isfield(psyinfo.range,'x') && numel(psyinfo.range.x > 2)
        nx = psyinfo.range.x(3);
    end

    % Prepare ranges
    if isfield(psyinfo.range,'mu')
        psy.mu(:,1,1) = linspace(psyinfo.range.mu(1),psyinfo.range.mu(2),nmu);
    else
        error('Cannot find a field for MU in initialization range struct.');
    end
    if isfield(psyinfo.range,'sigma')
        psy.logsigma(1,:,1) = linspace(log(psyinfo.range.sigma(1)), log(psyinfo.range.sigma(2)), nsigma);
    elseif isfield(psyinfo.range,'logsigma')
        psy.logsigma(1,:,1) = linspace(psyinfo.range.logsigma(1), psyinfo.range.logsigma(2), nsigma);
    else
        error('Cannot find a field for SIGMA in initialization range struct.');
    end
    if isfield(psyinfo.range,'lambda')
        psy.lambda(1,1,:) = linspace(psyinfo.range.lambda(1),psyinfo.range.lambda(2),nlambda);
    end
    if isfield(psyinfo,'x') && ~isempty(psyinfo.x)
        psy.x(1,1,1,:) = psyinfo.x(:);
    elseif isfield(psyinfo.range,'x') && ~isempty(psyinfo.range.x)
        psy.x(1,1,1,:) = linspace(psyinfo.range.x(1),psyinfo.range.x(2),nx);
    else
        error('Test grid X not provided in initialization struct.');
    end
end

% Default ranges
if isempty(psy.mu)
    psy.mu(:,1,1) = linspace(2,4,nmu);
end
if isempty(psy.logsigma)
    psy.logsigma(1,:,1) = linspace(log(0.01), log(1), nsigma);
end
if isempty(psy.lambda)
    psy.lambda(1,1,:) = linspace(0, 0.2, nlambda);
end
if isempty(psy.x)
    psy.x(1,1,1,:) = linspace(psy.mu(1),psy.mu(end),nx);
end

if isfield(psyinfo,'units')
    psy.units = psyinfo.units;
else
    psy.units.x = [];
    psy.units.mu = [];
    psy.units.sigma = []; 
    psy.units.lambda = []; 
end

% By default, wide Student's t prior on mu with slight preference for
% the middle of the stimulus range
muprior = [mean(psy.mu),0.5*(psy.mu(end)-psy.mu(1)),3];    % mean, sigma and nu
if isfield(psyinfo,'priors') && ~isempty(psyinfo.priors)
    if isfield(psyinfo.priors,'mu') && ~isempty(psyinfo.priors.mu)
        muprior(1:numel(psyinfo.priors.mu)) = psyinfo.priors.mu;        
    end
end
priormu = exp(logtpdf(psy.mu,muprior(1),muprior(2),muprior(3)));

% By default flat prior on log sigma (Jeffrey's 1/sigma prior in sigma 
% space); more in general log-Student-t prior
logsigmaprior = [mean(psy.logsigma),Inf,3];    % mean, sigma and nu
if isfield(psyinfo,'priors') && ~isempty(psyinfo.priors)
    if isfield(psyinfo.priors,'logsigma') && ~isempty(psyinfo.priors.logsigma)
        logsigmaprior(1:numel(psyinfo.priors.logsigma)) = psyinfo.priors.logsigma;
    end
end
priorlogsigma = exp(logtpdf(psy.logsigma,logsigmaprior(1),logsigmaprior(2),logsigmaprior(3)));

% Beta(a,b) prior on lambda, with correction
lambdaprior = [1,19];
if isfield(psyinfo,'priors') && ~isempty(psyinfo.priors)
    if isfield(psyinfo.priors,'lambda') && ~isempty(psyinfo.priors.lambda)
        lambdaprior = psyinfo.priors.lambda;        
    end
end

temp = psy.lambda(:)';
temp = [0, temp + 0.5*[diff(temp),0]];
a = lambdaprior(1); b = lambdaprior(2);
priorlambda(1,1,:) = betainc(temp(2:end),a,b) - betainc(temp(1:end-1),a,b);

% If using multiple conditions, divide prior
if Ncnd > 1
    fprintf('Sharing prior over LAMBDA across %d conditions.\n', Ncnd);
    priorlambda = priorlambda.^(1/Ncnd);
end

priormu = priormu./sum(priormu);
priorlogsigma = priorlogsigma./sum(priorlogsigma);
priorlambda = priorlambda./sum(priorlambda);

% Prior (posterior at iteration zero) over parameters
psy.post{1} = bsxfun(@times,bsxfun(@times,priormu,priorlogsigma),priorlambda);
for k = 2:Nfuns; psy.post{k} = psy.post{1}; end
for k = 1:Nfuns; psy.logupost{k} = log(psy.post{k}); end

% Posterior over lambda at iteration zero
psy.postlambda{1} = priorlambda;
for k = 2:Nfuns; psy.postlambda{k} = psy.postlambda{1}; end

% Define sigma in addition to log sigma
psy.sigma = exp(psy.logsigma);

% Guess rate for PCORRECT psychometric functions
if isfield(psyinfo,'gamma')
    psy.gamma = psyinfo.gamma;
else
    psy.gamma = [];
end

psy.f = [];

psy.xnext = [];

% Stimulus refractory time (avoid representing same stimulus for a while)
psy.reftime = []; psy.refradius = [];

% Refractory time constant (average run length of geometric distribution)
if isfield(psyinfo,'reftime')
    psy.reftime = psyinfo.reftime;
end
if isempty(psy.reftime); psy.reftime = 0; end

% Waiting radius (avoid representing stimuli within this radius)
if isfield(psyinfo,'refradius')
    psy.refradius = psyinfo.refradius;
end
if isempty(psy.refradius); psy.refradius = 0; end

% Initialize refractory times list
psy.reflist = zeros(size(psy.x));

% Set psychometric function
if isfield(psyinfo,'psychofun'); psy.psychofun = psyinfo.psychofun; end
[psy,Nfuns] = psyfunset(psy);

% Prior over psychometric functions
if Nfuns > 1
    if isfield(psyinfo,'psychoprior')
        psy.psychoprior = psyinfo.psychoprior;
    else
        psy.psychoprior = [];
    end    
    % Uniform prior by default
    if isempty(psy.psychoprior)
        psy.psychoprior = ones(1,Nfuns)/Nfuns;
    end
end

end

%--------------------------------------------------------------------------
function y = logtpdf(x,mu,sigma,nu)
%LOGTPDF Log pdf of Student's t distribution.

if sigma == Inf     % Flat log prior
    y = zeros(size(x));
elseif nu == Inf    % Student's t with infinite degrees of freedom is Gaussian
    y = -0.5*log(2*pi*sigma^2) -0.5*((x-mu)/sigma).^2;
else
    y = gammaln(0.5*(nu+1)) - gammaln(0.5*nu) - 0.5*log(pi*nu*sigma^2) ...
        - 0.5*(nu+1) * (1 + 1/nu * ((x-mu)/sigma).^2);
end

end

