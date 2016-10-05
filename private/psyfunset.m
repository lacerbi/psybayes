function [psy,Nfuns] = psyfunset(psy)
%PSYFUNSET Set psychometric function.

if ~isfield(psy,'psychofun'); psy.psychofun = []; end
if iscell(psy.psychofun)
    Nfuns = numel(psy.psychofun);
    cellflag = 1;
else
    Nfuns = 1;
    cellflag = 0;
end

% Select psychometric function
if ~isempty(psy.psychofun)
    if cellflag
        for k = 1:Nfuns
            if ischar(psy.psychofun{k})
                % Do nothing
            elseif isa(psy.psychofun{k},'function_handle')
                psy.psychofun{k} = func2str(psy.psychofun{k});
            else
                error('PSY.psychofun should be a psychometric function (function handle or function name), or a cell array of psychometric funtions.');
            end
            
        end
    else
        if ischar(psy.psychofun)
            % Do nothing
        elseif isa(psy.psychofun,'function_handle')
            psy.psychofun = func2str(psy.psychofun);
        else
            error('PSY.psychofun should be a psychometric function (function handle or function name), or a cell array of psychometric funtions.');
        end
    end        
else
    % Default psychometric functions (for YES/NO or PCORRECT)
    if ~isempty(psy.gamma)
        psy.psychofun = '@(x_,mu_,sigma_,lambda_,gamma_) psyfun_pcorrect(x_,mu_,sigma_,lambda_,gamma_)';
    else
        psy.psychofun = '@(x_,mu_,sigma_,lambda_,gamma_) psyfun_yesno(x_,mu_,sigma_,lambda_)';
    end    
end

psy.nfuns = Nfuns;