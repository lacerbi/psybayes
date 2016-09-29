function psybayes_plot(psy,refparams)
%PSYBAYES_PLOT Plot psychometric function and posterior from PSYBAYES.
%  PSYBAYES_PLOT(PSY) plots psychometric function and posterior over the
%  parameters as returned by PSYBAYES.
%
%  PSYBAYES_PLOT(PSY,REFPARAMS) also plots the reference parameter values
%  REFPARAMS=[MU,SIGMA,LAMBDA].
%
%  See also PSYBAYES, PSYTEST.

if nargin < 2; refparams = []; end

% Get psychometric function
psychofun = str2func(psy.psychofun);

% Arrange figure panels in a 2 x 2 vignette
rows = 2; cols = 2;
% rows = 1; cols = 4;   % Alternative horizontal arrangement

x = psy.x(:)';
xnext = psy.xnext;

% Plot psychometric function
subplot(rows,cols,1);

psimean = zeros(1,numel(x));    psisd = zeros(1,numel(x));
post = psy.post(:);
for ix = 1:numel(x)
    f = psychofun(x(ix),psy.mu,psy.sigma,psy.lambda,psy.gamma);
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
if ~isempty(refparams)
    psitrue = zeros(1,numel(x));
    for ix = 1:numel(x)
        psitrue(ix) = psychofun(x(ix),refparams(1),refparams(2),refparams(3),psy.gamma);
    end
    plot(x, psitrue, 'k','LineWidth',2);
end


if ~isempty(xnext)
    plot([xnext,xnext],[0,1],':r', 'LineWidth', 2);
end
if ~isempty(psy.data)
    scatter(psy.data(:,1),psy.data(:,2),20,'ko','MarkerFaceColor','r','MarkerEdgeColor','none');
end
box off; set(gca,'TickDir','out');    
if ~isempty(psy.units.x); string = [' (' psy.units.x ')']; else string = []; end
xlabel(['x' string]);
if isempty(psy.gamma)
    ylabel('Pr(response = 1)');
else
    ylabel('Pr(response correct)');        
end
axis([min(x) max(x), 0 1]);
title(['Psychometric function (trial ' num2str(num2str(psy.ntrial)) ')']);

% Plot posterior for mu
if numel(psy.mu) > 1
    subplot(rows,cols,3);
    y = sum(sum(psy.post,2),3);
    y = y/sum(y*diff(psy.mu(1:2)));    
    hold off;
    plot(psy.mu(:), y(:), 'k', 'LineWidth', 1);
    hold on;
    box off; set(gca,'TickDir','out');
    if ~isempty(psy.units.mu); string = [' (' psy.units.mu ')']; else string = []; end
    xlabel(['\mu' string]);
    ylabel('Posterior probability');
    % Compute SD of the posterior
    y = y/sum(y);
    ymean = sum(y.*psy.mu);
    ysd = sqrt(sum(y.*psy.mu.^2) - ymean^2);    
    if ~isempty(psy.units.mu); string = [' ' psy.units.mu]; else string = []; end
    title(['Posterior \mu = ' num2str(ymean,'%.2f') ' ± ' num2str(ysd,'%.2f') string])
    yl = get(gca,'Ylim'); axis([get(gca,'Xlim'),0,yl(2)]);
    if ~isempty(xnext)
        plot([xnext,xnext],get(gca,'Ylim'),':r', 'LineWidth', 2);
    end
end

% Plot posterior for sigma
if numel(psy.sigma) > 1
    subplot(rows,cols,2); hold off;
    y = sum(sum(psy.post,1),3);
    y = y/sum(y*diff(psy.sigma(1:2)));
    plot(psy.sigma(:), y(:), 'k', 'LineWidth', 1); hold on;
    box off; set(gca,'TickDir','out','XScale','log');
    %box off; set(gca,'TickDir','out','XScale','log','XTickLabel',{'0.1','1','10'});
    if ~isempty(psy.units.sigma); string = [' (' psy.units.sigma ')']; else string = []; end
    xlabel(['\sigma' string]);
    ylabel('Posterior probability');
    % title(['Marginal posterior distributions (trial ' num2str(tab.ntrial) ')']);
    % Compute SD of the posterior
    y = (y.*psy.sigma)/sum(y.*psy.sigma);
    ymean = sum(y.*psy.sigma);
    ysd = sqrt(sum(y.*psy.sigma.^2) - ymean^2);
    if ~isempty(psy.units.sigma); string = [' ' psy.units.sigma]; else string = []; end
    title(['Posterior \sigma = ' num2str(ymean,'%.2f') ' ± ' num2str(ysd,'%.2f') string]);
    yl = get(gca,'Ylim'); axis([psy.sigma(1),psy.sigma(end),0,yl(2)]);
end

% Plot posterior for lambda
if numel(psy.lambda) > 1
    subplot(rows,cols,4); hold off;
    y = sum(sum(psy.post,1),2);    
    y = y/sum(y*diff(psy.lambda(1:2)));
    plot(psy.lambda(:), y(:), 'k', 'LineWidth', 1); hold on;
    box off; set(gca,'TickDir','out');
    if ~isempty(psy.units.lambda); string = [' (' tabs.units.lambda ')']; else string = []; end
    xlabel(['\lambda' string]);
    ylabel('Posterior probability');    
    % Compute SD of the posterior
    y = y/sum(y);    
    ymean = sum(y.*psy.lambda);
    ysd = sqrt(sum(y.*psy.lambda.^2) - ymean^2);    
    if ~isempty(psy.units.lambda); string = [' ' tabs.units.lambda]; else string = []; end
    title(['Posterior \lambda = ' num2str(ymean,'%.2f') ' ± ' num2str(ysd,'%.2f') string])
    yl = get(gca,'Ylim'); axis([get(gca,'Xlim'),0,yl(2)]);
end

set(gcf,'Color','w');

% Add bars of reference parameterse
if ~isempty(refparams)
    if numel(psy.mu) > 1
        subplot(rows,cols,3);
        plot(refparams(1)*[1 1],get(gca,'Ylim'),'k','LineWidth',2);
        legendpanel = 3;
    end
    if numel(psy.lambda) > 1
        subplot(rows,cols,4);
        plot(refparams(3)*[1 1],get(gca,'Ylim'),'k','LineWidth',2);
        legendpanel = 4;
    end
    if numel(psy.sigma) > 1
        subplot(rows,cols,2);
        plot(refparams(2)*[1 1],get(gca,'Ylim'),'k','LineWidth',2);
        legendpanel = 2;
    end
    % Plot legend
    subplot(rows,cols,legendpanel);    
    axis([get(gca,'Xlim'),get(gca,'Ylim')]);
    hl = plot([-1 -1],[-1 -1],'k','LineWidth',2);
    h = legend(hl, 'Reference parameter value');
    set(h,'Location','NorthEast','Box','off');
end

end