function psy = psyfunset(psy)
%PSYFUNSET Set psychometric function.

% Select psychometric function
if isfield(psy,'psychofun') && ~isempty(psy.psychofun)
    if ischar(psy.psychofun)
        % Do nothing
    elseif isa(psy.psychofun,'function_handle')
        psy.psychofun = func2str(psy.psychofun);
    else
        error('PSY.psychofun should be a function handle or function name string for a psychometric function.');        
    end        
else
    % Default psychometric functions (for YES/NO or PCORRECT)
    if ~isempty(psy.gamma)
        psy.psychofun = '@(x_,mu_,sigma_,lambda_,gamma_) psyfun_pcorrect(x_,mu_,sigma_,lambda_,gamma_)';
    else
        psy.psychofun = '@(x_,mu_,sigma_,lambda_,gamma_) psyfun_yesno(x_,mu_,sigma_,lambda_)';
    end    
end