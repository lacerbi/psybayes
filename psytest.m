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
    % psy.gamma = [];   % Leave it empty for YES/NO psychometric functions

    % You can specify one or more user-defined psychometric functions (as a string)
    psy.psychofun{1} = '@(x,mu,sigma,lambda,gamma) psyfun_pcorrect(x,mu,sigma,lambda,gamma,@psynormcdf);';
    psy.psychofun{2} = '@(x,mu,sigma,lambda,gamma) psyfun_pcorrect(x,mu,sigma,lambda,gamma,@psygumbelcdf);';
    
    % Define range for stimulus and for parameters of the psychometric function
    % (lower bound, upper bound, number of points)
    psy.range.x = [1.5,4.5,61];
    psy.range.mu = [2,4,41];
    psy.range.sigma = [0.05,1,21];      % The range for sigma is automatically converted to log spacing
    psy.range.lambda = [0,0.1,15];
    %psy.range.lambda = [0.05,0.05,1];  % This would fix the lapse rate to 0.05

    % Define priors over parameters
    psy.priors.mu = [3,2];                  % mean and std of (truncated) Gaussian prior over MU
    psy.priors.logsigma = [log(0.5),Inf];   % mean and std of (truncated) Gaussian prior over log SIGMA (Inf std means flat prior)
    psy.priors.lambda = [1 50];             % alpha and beta parameter of beta pdf over LAMBDA

    % Units -- used just for plotting in axis labels and titles
    psy.units.x = 'cm';
    psy.units.mu = 'cm';
    psy.units.sigma = 'cm';
    psy.units.lambda = [];
    
    % Refractory time before presenting same stimulus again
    psy.reftime = 3;        % Expected number of trials (geometric distribution)
    psy.refradius = 0;      % Refractory radius around stimulus (in x units)
else    
    filename = 'psy.mat';   % Choose your file name
    load(filename,'psy');   % Load PSY structure from file
end

method = 'ent';     % Minimize the expected posterior entropy
% vars = [0 1 0];   % Minimize posterior entropy of the mean only
vars = [1 1 0];     % Minimize joint posterior entropy of mean and sigma (PSI-marginal method)
plotflag = 1;       % Plot visualization

%--------------------------------------------------------------------------

% Parameters of simulated observer
mu = 3.05;
sigma = 0.2;
lambda = 0.05;
gamma = psy.gamma;

%--------------------------------------------------------------------------
% Start running

string = ['Simulating an observer with MU = ' num2str(mu)];
if ~isempty(psy.units.mu); string = [string ' ' psy.units.mu]; end
string = [string '; SIGMA = ' num2str(sigma)];
if ~isempty(psy.units.sigma); string = [string ' ' psy.units.sigma]; end
string = [string '; LAMBDA = ' num2str(lambda)];
if ~isempty(psy.units.lambda); string = [string ' ' psy.units.lambda]; end
string = [string '.'];
display(string);
display('Press a key to simulate a trial.')

% Get first recommended point under the chosen optimization method
% (last argument is a flag for plotting)
[x,psy] = psybayes(psy, method, vars, [], []);

% Get psychometric function (only for simulating observer's responses)
psychofun = str2func(psy.psychofun{1});

pause;

Ntrials = 500;
EasyRate = 0.1;    % Frequency of easy trials

for iTrial = 1:Ntrials
    
    if rand < EasyRate
        % Simulate observer's response to an easy trial (only lapse/guessing)
        x = Inf;
        if isempty(gamma)
            r = rand < 1-lambda/2;
        else
            r = rand < 1-lambda*(1-gamma);            
        end
    else
        % Simulate observer's response given stimulus x
        r = rand < psychofun(x,mu,sigma,lambda,gamma);        
    end

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

