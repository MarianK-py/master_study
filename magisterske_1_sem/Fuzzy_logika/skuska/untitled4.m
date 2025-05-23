clf

data = rand([100 2]);

options = fcmOptions(...
    NumClusters=3,...
    Exponent=4,...
    Verbose=false);
[centers,U] = fcm(data,options);

maxU = max(U);
index1 = find(U(1,:) == maxU);
index2 = find(U(2,:) == maxU);
index3 = find(U(3,:) == maxU);

subplot(1,3,1)
[xq,yq] = meshgrid(0:.05:1, 0:.05:1);
vq = griddata(data(:,1),data(:,2),U(1,:),xq,yq);
surf(vq)
title("Cluster 1")

subplot(1,3,2)
[xq,yq] = meshgrid(0:.05:1, 0:.05:1);
vq = griddata(data(:,1),data(:,2),U(2,:),xq,yq);
surf(vq)
title("Cluster 2")

subplot(1,3,3)
[xq,yq] = meshgrid(0:.05:1, 0:.05:1);
vq = griddata(data(:,1),data(:,2),U(3,:),xq,yq);
surf(vq)
title("Cluster 3")