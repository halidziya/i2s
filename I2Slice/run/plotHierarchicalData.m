function plotHierarchicalData(X,slabels,labels,sampleindex)
    colormap(linspecer(max(max(labels))+1));
    if (~exist('sampleindex','var'))
        sampleindex = size(labels,1)
    end
    slabels = slabels(sampleindex,:)';
    labels = labels(sampleindex,:)';
    %labels = align_labels(labels');
    %slabels = align_labels(slabels');
    cla;
    scatter(X(:,1),X(:,2),10,slabels);
    for i=1:max(labels)
        if (sum(labels==i) > (size(X,2)+1)) 
        plot_gaussian_ellipsoid(mean(X(labels==i,1:2)),cov(X(labels==i,1:2)),'--',[0.5 0.5 0.5],2,0.5)
        end
    end

    for i=1:max(slabels)
        if (sum(slabels==i) > (size(X,2)+1))
        plot_gaussian_ellipsoid(mean(X(slabels==i,1:2)),cov(X(slabels==i,1:2)),'-',[0 0 0],2,2)
        end
    end

end