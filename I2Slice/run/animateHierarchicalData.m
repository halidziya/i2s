function animateHierarchicalData(X,slabels,labels)
filename = 'lastanim.gif';
for i=1:size(labels,1)
    plotHierarchicalData(X,slabels,labels,i)
    axis([min(X(:,1))-0.01 max(X(:,1))+0.01 min(X(:,2))-0.01 max(X(:,2))+0.01])
    mus=[];
    for j=1:max(labels(i,:))
        mus(j,:) = mean(X(labels(i,:)==j,:));
    end
    scatter(mus(:,1),mus(:,2),30,[0,0,0],'*');

    mus=[];
    for j=1:max(slabels(i,:))
        mus(j,:) = mean(X(slabels(i,:)==j,:));
    end
    scatter(mus(:,1),mus(:,2),50,[0,0,0],'o');
    scatter(mus(:,1),mus(:,2),30,[0,0,0]);
    if (size(mus,1)>2)
        voronoi(mus(:,1),mus(:,2))
    end
    drawnow limitrate
    frame = getframe(1);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
      if i == 1;
      imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
     else
      imwrite(imind,cm,filename,'gif','WriteMode','append');
      end
      
end
end