%PSYTEST Test and practical reference of PSYBAYES
%
%  Run PSYTEST for a demo of PSYBAYES.
%  Read the comments to learn the usage of PSYBAYES.
%
%  See also PSYBAYES, PSYBAYES_PLOT.

%--------------------------------------------------------------------------
% Definitions for PSYBAYES
newsession = 1;

if newsession
    % Initialize PSY structure
    psy = [];

    % Set chance level (for PCORRECT psychometric functions)
    psy.gamma = 0.5;        
    % psyinit.gamma = [];   % Leave it empty for YES/NO psychometric functions

    % Define range for stimulus and for parameters of the psychometric function
    % (lower bound, upper bound, number of points)
    psy.range.x = [1.5,4.5,61];
    psy.range.mu = [2,4,51];
    psy.range.sigma = [0.05,1,25];      % The range for sigma is automatically converted to log spacing
    psy.range.lambda = [0,0.4,25];

    % Define priors over parameters
    psy.priors.mu = [3,2];                  % mean and std of (truncated) Gaussian prior over MU
    psy.priors.logsigma = [log(0.5),Inf];   % mean and std of (truncated) Gaussian prior over log SIGMA (Inf std means flat prior)
    psy.priors.lambda = [1 19];             % alpha and beta parameter of beta pdf over LAMBDA

    % Units -- used just for plotting in axis labels and titles
    psy.units.x = 'cm';
    psy.units.mu = 'cm';
    psy.units.sigma = 'cm';
    psy.units.lambda = [];
else    
    filename = 'psy.mat';   % Choose your file name
    load(filename,'psy');   % Load PSY structure from file
end

method = 'ent';     % Minimize the expected posterior entropy
% vars = [1 0 0];   % Minimize posterior entropy of the mean only
vars = [1 1 1];     % Minimize joint posterior entropy of mean, sigma and lambda
plotflag = 1;       % Plot visualization

%--------------------------------------------------------------------------

% Parameters of simulated observer
mu = 3.05;
sigma = 0.2;
lambda = 0.05;

% Psychometric function for the simulated observer
if ~isfield(psy,'gamma') || isempty(psy.gamma)
    psychofun = @(x) lambda/2 + (1-lambda).*0.5*(1+erf((x-mu)./(sqrt(2)*sigma)));
else
    gamma = psy.gamma;
    psychofun = @(x) gamma + (1-gamma-lambda).*0.5*(1+erf((x-mu)./(sqrt(2)*sigma)));    
end

%--------------------------------------------------------------------------
% Start running

display(['Simulating an observer with MU = ' num2str(mu) ' cm; SIGMA = ' num2str(sigma) ' cm; LAMBDA = ' num2str(lambda) '.']);
display('Press a key to simulate a trial.')

% Get first recommended point under the chosen optimization method
% (last argument is a flag for plotting)
[x,psy] = psybayes(psy, method, vars, [], []);

pause;

Ntrials = 500;

for iTrial = 1:Ntrials
    % Simulate observer's response given stimulus x
    r = rand < psychofun(x);

    % Get next recommended point that minimizes predicted entropy 
    % given the current posterior and response r at x
    tic
    [x, psy] = psybayes(psy, method, vars, x, r);
    toc

    if plotflag
        trueparams = [mu,sigma,lambda];
        psybayes_plot(psy,trueparams);
    end
    
    drawnow;    
    % pause;
end

% Once you are done, clean PSY struct from temporary variables
[~,psy] = psybayes(psy);

% Save PSY struct, can be used in following sessions
% save(filename,'psy');

