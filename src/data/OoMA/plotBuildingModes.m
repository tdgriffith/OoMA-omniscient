function plotBuildingModes(Phi)
clf
for i = 1:size(Phi,2)
    subplot(1,size(Phi,2),i)
    z = Phi(:,i);
    phase = [cos(angle(z)) sin(angle(z))]; % the angles of the mode shape vector
    % choose the mode shape phase that maximizes the displacement
    [~,idx] = max([sum(abs(phase(:,1))) sum(abs(phase(:,2)))]);
    shape = abs(z).*phase(:,idx);
    plot([0; shape],0:size(Phi,1),'.-','markersize',15)
    hold on
    % plot the 'other side' of the mode shape
    plot([0; -shape],0:size(Phi,1),'.-','markersize',15)
    % plot the building centerline
    line([0 0],[0 size(Phi,1)],'color','k','linestyle','--')
    hold off
    yticks(1:size(Phi,1))
    xticks([])
    % remove axis ticks and resize axes
    inset =  get(gca,'tightinset');
    pos = get(gca,'position');
    set(gca,'position',[pos(1)-inset(1) pos(2)-inset(2) pos(3)+inset(1) pos(4)+inset(2)])
    if i == 1
        ylabel('Floor No.') % only label the x-axis at the first plot
    else
        yticklabels([])
    end
    xlabel(['Mode ' num2str(i)])
    axis tight
end
end