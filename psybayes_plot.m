function psybayes_plot(psystruct,refparams)
%PSYBAYES_PLOT Plot psychometric function and posterior from PSYBAYES.
%  PSYBAYES_PLOT(PSY) plots psychometric function and posterior over the
%  parameters as returned by PSYBAYES.
%
%  PSYBAYES_PLOT(PSY,REFPARAMS) also plots the reference parameter values
%  REFPARAMS=[MU,SIGMA,LAMBDA].
%
%  See also PSYBAYES, PSYTEST.

if nargin < 2; refparams = []; end

% Choose correct psychometric function (YES/NO or PCORRECT)
if ~isempty(psystruct.gamma)
    psychofun = @(x_,mu_,sigma_,lambda_) psyfun_pcorrect(x_,mu_,sigma_,lambda_,psystruct.gamma);
else    
    psychofun = @(x_,mu_,sigma_,lambda_) psyfun_yesno(x_,mu_,sigma_,lambda_);
end

% Arrange figure panels in a 2 x 2 vignette
rows = 2; cols = 2;
% rows = 1; cols = 4;   % Alternative horizontal arrangement

x = psystruct.x(:)';
xnext = psystruct.xnext;

% Plot psychometric function
subplot(rows,cols,1);

psimean = zeros(1,numel(x));    psisd = zeros(1,numel(x));
post = psystruct.post(:);
for ix = 1:numel(x)
    f = psychofun(x(ix),psystruct.mu,psystruct.sigma,psystruct.lambda);
    psimean(ix) = sum(f(:).*post);
    psisd(ix) = sqrt(sum(f(:).^2.*post) - psimean(ix)^2);
end    
hold off;
%area(x, psimean + psisd, 'EdgeColor', 'none', 'FaceColor', 0.8*[1 1 1]);
%hold on;
%area(x, psimean - psisd, 'EdgeColor', 'none', 'FaceColor', [1 1 1]);
fill([x fliplr(x)], [psimean+psisd, fliplr(psimean-psisd)], 0.8*[1 1 1], 'EdgeColor', 'none'); %, 'FaceColor', 0.8*[1 1 1]);
hold on;
plot(x, psimean,'k','LineWidth',1);

if ~isempty(xnext)
    plot([xnext,xnext],[0,1],':r', 'LineWidth', 2);
end
if ~isempty(psystruct.data)
    scatter(psystruct.data(:,1),psystruct.data(:,2),20,'ko','MarkerFaceColor','r','MarkerEdgeColor','none');
end
box off; set(gca,'TickDir','out');    
if ~isempty(psystruct.units.x); string = [' (' psystruct.units.x ')']; else string = []; end
xlabel(['x' string]);
if isempty(psystruct.gamma)
    ylabel('Pr(response = 1)');
else
    ylabel('Pr(response correct)');        
end
axis([min(x) max(x), 0 1]);
title(['Psychometric function (trial ' num2str(num2str(psystruct.ntrial)) ')']);

% Plot posterior for mu
subplot(rows,cols,3);
y = sum(sum(psystruct.post,2),3);
y = y/sum(y*diff(psystruct.mu(1:2)));    
hold off;
plot(psystruct.mu(:), y(:), 'k', 'LineWidth', 1);
hold on;
box off; set(gca,'TickDir','out');
if ~isempty(psystruct.units.mu); string = [' (' psystruct.units.mu ')']; else string = []; end
xlabel(['\mu' string]);
ylabel('Posterior probability');
% Compute SD of the posterior
y = y/sum(y);
ymean = sum(y.*psystruct.mu);
ysd = sqrt(sum(y.*psystruct.mu.^2) - ymean^2);    
if ~isempty(psystruct.units.mu); string = [' ' psystruct.units.mu]; else string = []; end
title(['Posterior \mu = ' num2str(ymean,'%.2f') ' ± ' num2str(ysd,'%.2f') string])
yl = get(gca,'Ylim'); axis([get(gca,'Xlim'),0,yl(2)]);
if ~isempty(xnext)
    plot([xnext,xnext],get(gca,'Ylim'),':r', 'LineWidth', 2);
end

% Plot posterior for sigma
subplot(rows,cols,2); hold off;
y = sum(sum(psystruct.post,1),3);
y = y/sum(y*diff(psystruct.sigma(1:2)));
plot(psystruct.sigma(:), y(:), 'k', 'LineWidth', 1); hold on;
box off; set(gca,'TickDir','out','XScale','log');
%box off; set(gca,'TickDir','out','XScale','log','XTickLabel',{'0.1','1','10'});
if ~isempty(psystruct.units.sigma); string = [' (' psystruct.units.sigma ')']; else string = []; end
xlabel(['\sigma' string]);
ylabel('Posterior probability');
% title(['Marginal posterior distributions (trial ' num2str(tab.ntrial) ')']);
% Compute SD of the posterior
y = (y.*psystruct.sigma)/sum(y.*psystruct.sigma);
ymean = sum(y.*psystruct.sigma);
ysd = sqrt(sum(y.*psystruct.sigma.^2) - ymean^2);
if ~isempty(psystruct.units.sigma); string = [' ' psystruct.units.sigma]; else string = []; end
title(['Posterior \sigma = ' num2str(ymean,'%.2f') ' ± ' num2str(ysd,'%.2f') string]);
yl = get(gca,'Ylim'); axis([psystruct.sigma(1),psystruct.sigma(end),0,yl(2)]);

% Plot posterior for lambda
subplot(rows,cols,4); hold off;
y = sum(sum(psystruct.post,1),2);    
y = y/sum(y*diff(psystruct.lambda(1:2)));
plot(psystruct.lambda(:), y(:), 'k', 'LineWidth', 1); hold on;
box off; set(gca,'TickDir','out');
if ~isempty(psystruct.units.lambda); string = [' (' tabs.units.lambda ')']; else string = []; end
xlabel(['\lambda' string]);
ylabel('Posterior probability');    
% Compute SD of the posterior
y = y/sum(y);    
ymean = sum(y.*psystruct.lambda);
ysd = sqrt(sum(y.*psystruct.lambda.^2) - ymean^2);    
if ~isempty(psystruct.units.lambda); string = [' ' tabs.units.lambda]; else string = []; end
title(['Posterior \lambda = ' num2str(ymean,'%.2f') ' ± ' num2str(ysd,'%.2f') string])
yl = get(gca,'Ylim'); axis([get(gca,'Xlim'),0,yl(2)]);

set(gcf,'Color','w');

% Add bars of reference parameterse
if ~isempty(refparams)
    subplot(rows,cols,3);
    plot(refparams(1)*[1 1],get(gca,'Ylim'),'k','LineWidth',2);
    subplot(rows,cols,2);
    plot(refparams(2)*[1 1],get(gca,'Ylim'),'k','LineWidth',2);
    axis([get(gca,'Xlim'),get(gca,'Ylim')]);
    hl = plot([-1 -1],[-1 -1],'k','LineWidth',2);
    h = legend(hl, 'True parameter value');
    set(h,'Location','NorthEast','Box','off');
    subplot(rows,cols,4);
    plot(refparams(3)*[1 1],get(gca,'Ylim'),'k','LineWidth',2);
end

end