clf

data = rand([100 2]);

options = fcmOptions(...
    NumClusters=3,...
    Exponent=4,...
    Verbose=false);
[centers,U] = fcm(data,options);

Thresh = 0.3;
maxU = max(U);
index1 = find(U(1,:) >= Thresh);
nonIndex1 = setdiff(1:100, index1);
index2 = find(U(2,:) >= Thresh);
nonIndex2 = setdiff(1:100, index2);
index3 = find(U(3,:) >= Thresh);
nonIndex3 = setdiff(1:100, index3);

subplot(1,3,1)
plot(data(index1,1),data(index1,2),'ob')
hold on
plot(data(nonIndex1,1),data(nonIndex1,2),'or')
hold on
plot(centers(1,1),centers(1,2),'xblack')

subplot(1,3,2)
plot(data(index2,1),data(index2,2),'ob')
hold on
plot(data(nonIndex2,1),data(nonIndex2,2),'or')
hold on
plot(centers(2,1),centers(2,2),'xblack')

subplot(1,3,3)
plot(data(index3,1),data(index3,2),'ob')
hold on
plot(data(nonIndex3,1),data(nonIndex3,2),'or')
hold on
plot(centers(3,1),centers(3,2),'xblack')