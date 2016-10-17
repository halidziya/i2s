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
    covs=[];
    for i=1:max(labels)
        if (sum(labels==i) > (size(X,2)+1)) 
            covs(i,:,:) = cov(X(labels==i,1:2));
            plot_gaussian_ellipsoid(mean(X(labels==i,1:2)),squeeze(covs(i,:,:)),'--',[0.4 0.4 0.4],2,1)
        else
            covs(i,:,:) = zeros(2);
        end
    end

    for i=1:max(slabels)
        if (sum(slabels==i) > (size(X,2)+1))
        plot_gaussian_ellipsoid(mean(X(slabels==i,1:2)),cov(X(slabels==i,1:2)),'-',[0 0 0],2,2)
        end
    end
    lmap = unique([slabels labels],'rows');
    for i=1:max(slabels)
        if (sum(slabels==i) > (size(X,2)+1))
            try
                plot_gaussian_ellipsoid(mean(X(slabels==i,1:2)),squeeze(mean(covs(lmap(lmap(:,1)==i,2),:,:),1)),'-',[0.4 0.2 0.2],2,2)
            catch 
                fprint('Skipped')
            end
        end
    end
end