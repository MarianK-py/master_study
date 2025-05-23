 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cluster Quasi-Random Data Using Fuzzy C-Means Clustering
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load fcmdata0.dat
plot(fcmdata(:,1),fcmdata(:,2),'o')
[center,U,objFcn] = fcm(fcmdata,2);

figure
plot(objFcn)
title('Objective Function Values')
xlabel('Iteration Count')
ylabel('Objective Function Value')

maxU = max(U);
index1 = find(U(1,:) == maxU);
index2 = find(U(2,:) == maxU);
figure
line(fcmdata(index1,1), fcmdata(index1,2), 'linestyle','none','marker','o','color','g')
line(fcmdata(index2,1),fcmdata(index2,2),'linestyle','none','marker','x','color','r')
hold on
plot(center(1,1),center(1,2),'ko','markersize',15,'LineWidth',2)
plot(center(2,1),center(2,2),'kx','markersize',15,'LineWidth',2)    

figure
plot(U(1,:),'-b')
figure
plot(U(2,:),'-r')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adjust Fuzzy Overlap in Fuzzy C-Means Clustering
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rng('default');
data = rand(100,2);
M = [1.1 2.0 3.0 4.0];
for i = 1:4
% Cluster the data.
options = [M(i) NaN NaN 0];
[centers,U] = fcm(data,2,options);
% Classify the data points.
maxU = max(U);
index1 = find(U(1,:) == maxU);
index2 = find(U(2,:) == maxU);
% Find data points with lower maximum membership values.
index3 = find(maxU < 0.6);
% Calculate the average maximum membership value.
averageMax = mean(maxU);
% Plot the results.
subplot(2,2,i)
plot(data(index1,1),data(index1,2),'ob')
hold on
plot(data(index2,1),data(index2,2),'or')
plot(data(index3,1),data(index3,2),'xk','LineWidth',2)
plot(centers(1,1),centers(1,2),'xb','MarkerSize',15,'LineWidth',3)
plot(centers(2,1),centers(2,2),'xr','MarkerSize',15,'LineWidth',3)
hold off
title(["M = " num2str(M(i)) ", Ave. Max. = " num2str(averageMax,3)])
end
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fuzzy C-Means Clustering
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
load fcmdata3.dat
% writematrix(fcmdata3,'fcmdata3.dat')
dataset=fcmdata3;
N = 4;
exp = 2;
maxIter = 100;
minImprove = 0.00001;
displayObjective = false;
options = [exp maxIter minImprove displayObjective];
[C,U] = fcm(dataset,N,options);
 
maxU = max(U);
index = cell(N,1);
for i=1:N
index{i} = find(U(i,:) == maxU);
end
 
figure
hold on
for i=1:N
plot(dataset(index{i},1),dataset(index{i},2),'o')
plot(C(i,1),C(i,2),'xk','MarkerSize',15,'LineWidth',3)
end
hold off

cluster = 2;
[X,Y] = meshgrid(0:0.05:1, 0:0.05:1);
Z = griddata(dataset(:,1),dataset(:,2),U(cluster,:),X,Y);
surf(X,Y,Z)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fuzzy C-Means Clustering for Iris Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
load iris.dat
setosaIndex = iris(:,5)==1;
versicolorIndex = iris(:,5)==2;
virginicaIndex = iris(:,5)==3;
setosa = iris(setosaIndex,:);
versicolor = iris(versicolorIndex,:);
virginica = iris(virginicaIndex,:);

Characteristics = {'sepal length','sepal width','petal length','petal width'};
pairs = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
for i = 1:6
x = pairs(i,1);
y = pairs(i,2);
subplot(2,3,i)
plot([setosa(:,x) versicolor(:,x) virginica(:,x)],...
     [setosa(:,y) versicolor(:,y) virginica(:,y)], '.')
xlabel(Characteristics{x})
ylabel(Characteristics{y})
end

Nc = 3;
M = 2.0;
maxIter = 100;
minImprove = 1e-6;
clusteringOptions = [M maxIter minImprove true];
[centers,U] = fcm(iris,Nc,clusteringOptions);

for i = 1:6
subplot(2,3,i)
 for j = 1:Nc
 x = pairs(i,1);
 y = pairs(i,2);
 text(centers(j,x),centers(j,y),int2str(j),'FontWeight','bold');
 end
end
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model Suburban Commuting Using Subtractive Clustering and ANFIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load trafficData.mat
subplot(2,1,1)
plot(datain)
legend('population','number of dwellings','vehicle ownership',...
       'median income','total employment','Location','northwest')
title('Input Variables')
xlabel('Data Point')
subplot(2,1,2)
plot(dataout)
legend('number of trips')
title('Output Variable')
xlabel('Data Point')

[C,S] = subclust([datain dataout],0.5);

figure
plot(datain(:,5),dataout(:,1),'.',C(:,5),C(:,6),'r*')
legend('Data points','Cluster centers','Location','southeast')
xlabel('Total Employment')
ylabel('Number of Trips')
title('Data and Cluster Centers')

% S = 1×6
% 1.1621 0.4117 0.6555 7.6139 2.8931 1.4395

fisOpt = genfisOptions('SubtractiveClustering','ClusterInfluenceRange',0.5);
fis = genfis(datain,dataout,fisOpt);

showrule(fis,'Format','symbolic')
figure
plotmf(fis,'input',1)

fuzout = evalfis(fis,datain);
trnRMSE = norm(fuzout-dataout)/sqrt(length(fuzout));
valfuzout = evalfis(fis,valdatain);
valRMSE = norm(valfuzout-valdataout)/sqrt(length(valfuzout));

figure
plot(valdataout)
hold on
plot(valfuzout,'o')
hold off
ylabel('Output value')
legend('Validation data','FIS output','Location','northwest')

opt = anfisOptions('InitialFIS',fis,'EpochNumber',20,'InitialStepSize',0.1);
fis2 = anfis([datain dataout],opt);
 
fuzout2 = evalfis(fis2,datain);
trnRMSE2 = norm(fuzout2-dataout)/sqrt(length(fuzout2));
valfuzout2 = evalfis(fis2,valdatain);
valRMSE2 = norm(valfuzout2-valdataout)/sqrt(length(valfuzout2));

figure
plot(valdataout)
hold on
plot(valfuzout,'o')
plot(valfuzout2,'x')
hold off
ylabel('Output value')
legend('Validation data',...
       "Initial FIS: RMSE = " + num2str(valRMSE), ...
       "Tuned FIS: RMSE = " + num2str(valRMSE2), ...
       'Location','northwest')

opt.EpochNumber = 200;
opt.ValidationData = [valdatain valdataout];
opt.DisplayANFISInformation = 0;
opt.DisplayErrorValues = 0;
opt.DisplayStepSize = 0;
opt.DisplayFinalResults = 0;
[fis3,trnErr,stepSize,fis4,valErr] = anfis([datain dataout],opt);

% fis3 is the FIS object when the training error reaches a minimum.
% fis4 is the snapshot FIS object when the validation data error reaches a minimum.
% stepSize is a history of the training step sizes.
% trnErr is the RMSE using the training data.
% valErr is the RMSE using the validation data for each training epoch.

fuzout4 = evalfis(fis4,datain);
trnRMSE4 = norm(fuzout4-dataout)/sqrt(length(fuzout4));
valfuzout4 = evalfis(fis4,valdatain);
valRMSE4 = norm(valfuzout4-valdataout)/sqrt(length(valfuzout4));

figure
plot(valdataout)
hold on
plot(valfuzout2,'o')
plot(valfuzout4,'x')
hold off
ylabel('Output value')
legend('Validation data',...
       "Tuned FIS: RMSE = " + num2str(valRMSE2), ...
       "Min val. error FIS: RMSE = " + num2str(valRMSE4), ...
       'Location','northwest')

figure
plot(trnErr)
title('Training Error')
xlabel('Number of Epochs')
ylabel('Error')

figure
plot(valErr)
title('Validation Error')
xlabel('Number of Epochs')
ylabel('Error')

findcluster
load clusterdemo.dat


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HELP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rng('default')
xy = -2.5 + 5*rand([200 2]);
x = xy(:,1);
y = xy(:,2);
v = x.*exp(-x.^2-y.^2);
[xq,yq] = meshgrid(-2:.2:2, -2:.2:2);
vq = griddata(x,y,v,xq,yq);


options = [exp maxIter minImprove displayObjective];
[centers,U] = fcm(data,Nc,options);
opt = genfisOptions('FCMClustering');
opt.NumClusters = Nc;
opt.Exponent = options(1);
opt.MaxNumIteration = options(2);
opt.MinImprovement = options(3);
opt.Verbose = options(4);
inputData = data(:,1:M);
outputData = data(:,M+1:end);
fis = genfis(inputData,outputData,opt);


subclust
Find cluster centers using subtractive clustering

centers = subclust(data,clusterInfluenceRange)
centers = subclust(data,clusterInfluenceRange,Name,Value)
[centers,sigma] = subclust( ___ )

centers = subclust(data,clusterInfluenceRange) 
clusters input data using subtractive clustering with the specified cluster 
influence range, and returns the computed cluster centers. 
The subtractive clustering algorithm estimates the number of clusters in the input data.
centers = subclust(data,clusterInfluenceRange,Name,Value) 
clusters data using algorithm options specified by one or more 
Name,Value arguments.
[centers,sigma] = subclust( ___ ) 
returns the sigma values specifying the range of influence of a cluster 
center in each of the data dimensions.

Find Cluster Centers Using Subtractive Clustering
load clusterDemo.dat
Find cluster centers using the same range of influence for all dimensions.
C = subclust(clusterDemo,0.6);
Each row of C contains one cluster center.
C
C = 3×3
0.5779 0.2355 0.5133
0.7797 0.8191 0.1801
0.1959 0.6228 0.8363

Specify Bounds for Subtractive Clustering
load clusterDemo.dat
Define minimum and maximum normalization bounds for each data dimension. 
Use the same bounds for each dimension.
dataScale = [-0.2 -0.2 -0.2;
              1.2 1.2 1.2];
Find cluster centers.
C = subclust(clusterDemo,0.5,'DataScale',dataScale);

Specify Options for Subtractive Clustering
load clusterDemo.dat
Specify the following clustering options:
Squash factor of 2.0 - Only find clusters that are far from each other.
Accept ratio 0.8 - Only accept data points with a strong potential for being cluster centers.
Reject ratio of 0.7 - Reject data points if they do not have a strong potential 
for being cluster centers.
Verbosity flag of 0 - Do not print progress information to the command window.
options = [2.0 0.8 0.7 0];
Find cluster centers, using a different range of influence for each 
dimension and the specified options.
C = subclust(clusterDemo,[0.5 0.25 0.3],'Options',options);

Obtain Cluster Influence Range for Each Data Dimension
load clusterDemo.dat
Cluster data, returning cluster sigma values, S.
[C,S] = subclust(clusterDemo,0.5);
Cluster sigma values indicate the range of influence of the computed 
cluster centers in each data dimension.


Input Arguments

data — Data set to be clustered M-by-N array
Data to be clustered, specified as an M-by-N array, where M is the number 
of data points and N is the number of data dimensions.

clusterInfluenceRange — Range of influence of the cluster center
scalar value in the range [0, 1] | vector
Range of influence of the cluster center for each input and output assuming 
the data falls within a unit hyperbox, specified as one of the following:
Scalar value in the range [0 1] — Use the same influence range for all inputs and outputs.
Vector — Use different influence ranges for each input and output.
Specifying a smaller range of influence usually creates more and smaller 
data clusters, producing more fuzzy rules.

Name-Value Pair Arguments
Specify optional pairs of arguments as Name1=Value1,...,NameN=ValueN, 
where Name is the argument name and Value is the corresponding value. 
Name-value arguments must appear after other arguments, but the order 
of the pairs does not matter.
Example: centers = subclust(data,0.5,DataScale=10)
Example: centers = subclust(data,0.5,"DataScale",10)

DataScale — Data scale factors
"auto" (default) | 2-by-N array
Data scale factors for normalizing input and output data into a unit hyperbox, 
specified as a 2-by-N array, where N is the total number of inputs and outputs. 
Each column of DataScale specifies the minimum value in the first row and 
the maximum value in the second row for the corresponding input or output data set.
When DataScale is "auto", the subclust function uses the actual minimum and maximum values
in the data to be clustered.

Options — Clustering options
vector
Clustering options, specified as a vector with the following elements.

Options(1) — Squash factor
1.25 (default) | positive scalar
Squash factor for scaling the range of influence of cluster centers, 
specified as a positive scalar. 
A smaller squash factor reduces the potential for outlying points to be 
considered as part of a cluster, which usually creates more and smaller data clusters.

Options(2) — Acceptance ratio
0.5 (default) | scalar value in the range [0, 1]
Acceptance ratio, defined as a fraction of the potential of the first 
cluster center, above which another data point is accepted as a cluster center, 
specified as a scalar value in the range [0, 1]. 
The acceptance ratio must be greater than the rejection ratio.

Options(3) — Rejection ratio
0.15 (default) | scalar value in the range [0, 1]
Rejection ratio, defined as a fraction of the potential of the first cluster center, 
below which another data point is rejected as a cluster center, 
specified as a scalar value in the range [0, 1]. 
The rejection ratio must be less than acceptance ratio.

Options(4) — Information display flag
false (default) | true
Information display flag indicating whether to display progress information 
during clustering, specified as one of the following:
false — Do not display progress information.
true — Display progress information.


Output Arguments

centers — Cluster centers J-by-N array
Cluster centers, returned as a J-by-N array, where J is the number of clusters and 
N is the number of data dimensions.

sigma — Range of influence of cluster centers N-element row vector
Range of influence of cluster centers for each data dimension, returned as 
an N-element row vector.
All cluster centers have the same set of sigma values.


To generate a fuzzy inference system using subtractive clustering, use 
the genfis command. 
For example, suppose you cluster your data using the following syntax:
C = subclust(data,clusterInfluenceRange,"DataScale",dataScale,"Options",options);
where the first M columns of data correspond to input variables, and the remaining columns
correspond to output variables.
You can generate a fuzzy system using the same training data and subtractive clustering
configuration. 
To do so:

Configure clustering options.
opt = genfisOptions("SubtractiveClustering");
opt.ClusterInfluenceRange = clusterInfluenceRange;
opt.DataScale = dataScale;
opt.SquashFactor = options(1);
opt.AcceptRatio = options(2);
opt.RejectRatio = options(3);
opt.Verbose = options(4);

Extract input and output variable data.
inputData = data(:,1:M);
outputData = data(:,M+1:end);

Generate FIS structure.
fis = genfis(inputData,outputData,opt);
The fuzzy system, fis, contains one fuzzy rule for each cluster, and each input and output
variable has one membership function per cluster. 
You can generate only Sugeno fuzzy systems using subtractive clustering. 


Algorithms

Subtractive clustering assumes that each data point is a potential cluster center. 
The algorithm does the following:
1 Calculate the likelihood that each data point would define a cluster center, based on the density
of surrounding data points.
2 Choose the data point with the highest potential to be the first cluster center.
3 Remove all data points near the first cluster center. The vicinity is determined using
clusterInfluenceRange.
4 Choose the remaining point with the highest potential as the next cluster center.
5 Repeat steps 3 and 4 until all the data is within the influence range of a cluster center.

