clf

data = rand([100 2]);

for i = [1 2 3 4; 1.1 2 4 8]
options = fcmOptions(...
    NumClusters=3,...
    Exponent=i(2),...
    Verbose=false);
[centers,U] = fcm(data,options);

maxU = max(U);
index1 = find(U(1,:) == maxU);
index2 = find(U(2,:) == maxU);
index3 = find(U(3,:) == maxU);

subplot(2,2,i(1))

plot(data(index1,1),data(index1,2),'ob')
hold on
plot(data(index2,1),data(index2,2),'or')
hold on
plot(data(index3,1),data(index3,2),'og')
hold on
plot(centers([1 2 3],1),centers([1 2 3],2),'xblack')
title("Exponent: "+i(2))
end