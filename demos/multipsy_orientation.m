%MULTIPSY_ORIENTATION Test of PSYBAYES for orientation experiment.
%
%  See also PSYBAYES, PSYBAYES_PLOT, PSYTEST.

%--------------------------------------------------------------------------
% Definitions for PSYBAYES
newsession = ~exist('psy','var') || isempty(psy);
Ncnd = 4;   % Different conditions

if newsession
    psy = [];
    for c = 1:Ncnd
        % Initialize PSY structure
        psy{c} = [];

        % Specify user-defined psychometric function (as a string)
        psy{c}.psychofun{1} = '@(x,mu,sigma,lambda,gamma) psyfun_yesno(x,mu,sigma,lambda,gamma,@psynormcdf);';

        % Define range for stimulus and for parameters of the psychometric function
        % (lower bound, upper bound, number of points)

        x = [0.25:0.25:1 1.4:0.4:3 3.5:0.5:5 6:10 12:2:20 22.5:2.5:45];
        x = [-fliplr(x),x];

        psy{c}.x = x;           % Stimulus grid in orientation degrees
        psy{c}.range.mu = [-12,12,31];     % Psychometric function mean
        psy{c}.range.sigma = [0.5,45,41];                   % The range for sigma is automatically converted to log spacing
        psy{c}.range.lambda = [0,0.5,21];
        % psy.range.lambda = [0.05-eps,0.05+eps,2];  % This would fix the lapse rate to 0.05

        % Define priors over parameters
        psy{c}.priors.mu = [0,2.5];              % mean and sigma of (truncated) Student's t prior over MU
        psy{c}.priors.logsigma = [2.05,1.40];  % mean and sigma of (truncated) Student's t prior over log SIGMA (Inf std means flat prior)
        psy{c}.priors.lambda = [1 19];         % alpha and beta parameters of beta pdf over LAMBDA

        % Units -- used just for plotting in axis labels and titles
        psy{c}.units.x = 'deg';
        psy{c}.units.mu = 'deg';
        psy{c}.units.sigma = 'log deg';
        psy{c}.units.lambda = [];
        psy{c}.units.psychofun = {'Normal'};

        % Refractory time before presenting same stimulus again
        psy{c}.reftime = 0;        % Expected number of trials (geometric distribution)
        psy{c}.refradius = 0;      % Refractory radius around stimulus (in x units)

        % Force stimulus placement to be left and right of mean psychometric curve
        psy{c}.forcesymmetry = 1;
    end
    
else    
    %filename = 'psy.mat';   % Choose your file name
    %load(filename,'psy');   % Load PSY structure from file
end

method = 'ent';     % Minimize the expected posterior entropy
vars = [1 1 1];     % Minimize posterior entropy of mean and sigma
% vars = [1 1 1];     % Minimize joint posterior entropy of mean, sigma and lambda
plotflag = 1;       % Plot visualization

%--------------------------------------------------------------------------

% Parameters of simulated observer
mu = [-1 1 0.5 2];
sigma = [exp(1.2) exp(1.5) exp(1.7) exp(2)];
lambda = [0.1 0.1 0.1 0.1];
gamma = [];

%--------------------------------------------------------------------------
% Start running

string = ['Simulating an observer with MU = ' num2str(mu)];
if ~isempty(psy{1}.units.mu); string = [string ' ' psy{1}.units.mu]; end
string = [string '; SIGMA = ' num2str(sigma)];
if ~isempty(psy{1}.units.sigma); string = [string ' ' psy{1}.units.sigma]; end
string = [string '; LAMBDA = ' num2str(lambda)];
if ~isempty(psy{1}.units.lambda); string = [string ' ' psy{1}.units.lambda]; end
string = [string '.'];
display(string);
display('Press a key to simulate a trial.')

% Get first recommended point under the chosen optimization method
% (last argument is a flag for plotting)
[x,psy] = psybayes_joint(psy, method, vars, [], [], 1);

% Get psychometric function (only for simulating observer's responses)
psychofun = cell(1,Ncnd);
for c = 1:Ncnd
    psychofun{c} = str2func(psy{1}.psychofun{1});
end

pause;

Ntrials = 500;
EasyRate = 0.0;    % Frequency of easy trials
c = 1;             % First condition

for iTrial = 1:Ntrials
    
    % Randomize condition after first trial
    if iTrial > 1; c = randi(Ncnd); end

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
        r = rand < psychofun{c}(x,mu(c),sigma(c),lambda(c),gamma);        
    end
    
    % Get next recommended point that minimizes predicted entropy 
    % given the current posterior and response r at x
    tic
    [x, psy] = psybayes_joint(psy, method, vars, x, r, c);
    toc

    if plotflag
        figure(c);
        trueparams = [mu(c),sigma(c),lambda(c)];
        psybayes_plot(psy{c},trueparams);
    end
    
    % Print posterior over distinct psychometric curves
    for k = 1:numel(psy{c}.psychofun)
        if k < numel(psy{c}.psychofun); sep = ', '; else sep = '.  '; end
        fprintf(['%s: %.3f' sep], psy{c}.units.psychofun{k}, psy{c}.psychopost(k));
    end
    drawnow;    
    % pause;
end

% Once you are done, clean PSY struct from temporary variables
[~,psy] = psybayes_joint(psy);

% Save PSY struct, can be used in following sessions
% save(filename,'psy');

