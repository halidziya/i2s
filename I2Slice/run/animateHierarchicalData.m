function animateHierarchicalData(X,slabels,labels)
for i=1:size(labels,1)
    plotHierarchicalData(X,slabels,labels,i)
    axis([min(X(:,1))-0.01 max(X(:,1))+0.01 min(X(:,2))-0.01 max(X(:,2))+0.01])
    drawnow limitrate
end
end