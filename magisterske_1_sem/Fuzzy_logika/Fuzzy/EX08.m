
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ANFIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fis = anfis(trainingData)
fis = anfis(trainingData,options)
[fis,trainError] = anfis( ___ )
[fis,trainError,stepSize] = anfis( ___ )
[fis,trainError,stepSize,chkFIS,chkError] = anfis(trainingData,options)

options — Training option.
Using options, you can specify:
An initial FIS structure to tune, options.InitialFIS.
Validation data for preventing overfitting to training data, options.ValidationData.
Training algorithm options, such as the maximum number of training epochs,
options.EpochNumber, or the training error goal, options.ErrorGoal.
Whether to display training progress information, such as the training error values for each
training epoch, options.DisplayErrorValues.

fis — Trained fuzzy inference system.
Trained fuzzy inference system with membership function parameters tuned using the training data,
returned as a mamfis or sugfis object. This fuzzy system corresponds to the epoch for which the
training error is smallest. If two epochs have the same minimum training error, the FIS from the
earlier epoch is returned.

chkFIS — Tuned FIS for which the validation error is minimum.
Tuned FIS for which the validation error is minimum, returned as a mamfis or sugfis object. If two
epochs have the same minimum validation error, the FIS from the earlier epoch is returned.
chkFIS is returned only when you specify validation data using options.ValidationData.

Generate and train a fuzzy inference system. 
By default, the FIS structure is created using a grid
partition of the input variable range with two membership functions.
load fuzex1trnData.mat
fis = anfis(fuzex1trnData);
%%
x = fuzex1trnData(:,1);
anfisOutput = evalfis(fis,x);
plot(x,fuzex1trnData(:,2),'*r',x,anfisOutput,'.b')
legend('Training Data','ANFIS Output','Location','NorthWest')
%%
opt = anfisOptions('InitialFIS',4,'EpochNumber',40);
fis = anfis(fuzex1trnData,opt);

x = (0:0.1:10)';
y = sin(2*x)./exp(x/5);

genOpt = genfisOptions('GridPartition');
genOpt.NumMembershipFunctions = 5;
genOpt.InputMembershipFunctionType = 'gaussmf';
inFIS = genfis(x,y,genOpt);

opt = anfisOptions('InitialFIS',inFIS);
opt.DisplayANFISInformation = 0;
opt.DisplayErrorValues = 0;
opt.DisplayStepSize = 0;
opt.DisplayFinalResults = 0;
outFIS = anfis([x y],opt);

plot(x,y,x,evalfis(outFIS,x))
legend('Training Data','ANFIS Output')

opt = anfisOptions('InitialFIS',4,'EpochNumber',40);
opt.ValidationData = fuzex2chkData;
[trainFIS,trainFISError,~,validationFIS,validationFISError] = anfis(fuzex2trnData,opt);
trainFISRMSE = min(trainFISError);
validationFISRMSE = min(validationFISError);

opt.StepSizeIncreaseRate = 2*opt.StepSizeIncreaseRate;
[fis,~,stepSize] = anfis([x y],opt);
plot(stepSize)

opt.ValidationData = fuzex1chkData;
[fis,trainError,stepSize,chkFIS,chkError] = anfis(fuzex1trnData,opt);
x = [1:40];
plot(x,trainError,'.b',x,chkError,'*r')
opt.ErrorGoal

x = (0:0.1:10)';
y = sin(2*x)./exp(x/5);
options = genfisOptions('GridPartition');
options.NumMembershipFunctions = 5;
fisin = genfis(x,y,options);
[in,out,rule] = getTunableSettings(fisin);
opt = tunefisOptions('Method','anfis');
fisout = tunefis(fisin,[in;out],x,y,opt);
%% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXERCISES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

%Uloha 1:
%Training:
x = (0:0.1:10)';
y = sin(2*x)./exp(x/5);

train1=[x y];

opt = anfisOptions('InitialFIS',4,'EpochNumber',60);
a1 = anfis(train1,opt);
plot(x,y,x,evalfis(a1,x))

%%
save train1.mat train1
load train1.mat train1
%%
Uloha 2:
Training:
x = (0:0.1:10)';
y = sin(2*x)./exp(x/5);
Validating:
x = (0.05:0.1:10.05)';
y = sin(2*x)./exp(x/5);

Uloha 3:
Training:
x = (0:0.05:2*pi-0.3)';
y = sin(x)./cos(x);
Validating:
x = (0.1:0.05:2*pi-0.2)';
y = sin(x)./cos(x);

%%
%Uloha 4:
%Training:
d = (0:0.05:2*pi);
[x y] = meshgrid(d,d); 
z = sin(x)+cos(y);
train4=[x(:) y(:) z(:)];
%opt = genfisOptions('GridPartition');
opt = genfisOptions("FCMClustering","NumClusters",10);
opt.MaxNumIteration = 1000;
%opt.NumMembershipFunctions = [10 10];
%opt.InputMembershipFunctionType = ["gaussmf" "gaussmf"];
ansugex4init = genfis([x(:) y(:)],[z(:)],opt);
ansugex4finOUT = evalfis(ansugex4init,[x(:) y(:)]);
ansugex4finGRIDOUT = griddata(x(:),y(:),ansugex4finOUT(:),x,y);
surf(x,y,ansugex4finGRIDOUT)
hold on
surf(x,y,z)
%%
writefis(ansugex4init,'ansugex4init.fis');
%%
Uloha 5:
Training:
d = (-10:0.1:10);
[x y] = meshgrid(d,d); 
r = (x.^2+y.^2)^(1/3);
z = sin(r)./r;
train5=[x(:) y(:) real(z(:))];
opt = genfisOptions('GridPartition');
opt.NumMembershipFunctions = [5 5];
opt.InputMembershipFunctionType = ["gaussmf" "gaussmf"];
ansugex5init = genfis([x(:) y(:)],[real(z(:))],opt);
writefis(ansugex5init,'ansugex5init.fis');

Uloha 6:
Training:
a = (0:0.05:2*pi);
ca = cos(a);
m = 0:0.05:10;
[x y] = meshgrid(m,ca); 
z = (x.^2).*sqrt(1-(y.^2));
% surf(x,y,z)
train6=[x(:) y(:) z(:)];
opt = genfisOptions('GridPartition');
opt.NumMembershipFunctions = [5 5];
opt.InputMembershipFunctionType = ["gaussmf" "gaussmf"];
ansugex6init = genfis([x(:) y(:)],[z(:)],opt);
writefis(ansugex6init,'ansugex6init.fis');
% ansugex6finOUT = evalfis(ansugex6fin,[x(:) y(:)])
% ansugex6finGRIDOUT = griddata(x(:),y(:),ansugex6finOUT(:),x,y)
% ansugex6finGRIDOUT = griddata(x(:),y(:),ansugex6finOUT,x,y)
% surf(x,y,ansugex6finGRIDOUT)
% colormap(cool) colormap(hot) alpha(.1)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TUNING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

evalfis
Evaluate fuzzy inference system

output = evalfis(fis,input)
output = evalfis(fis,input,options)
[output,fuzzifiedIn,ruleOut,aggregatedOut,ruleFiring] = evalfis( ___ )

output = evalfis(fis,input) 
evaluates the fuzzy inference system fis for the input values in
input and returns the resulting output values in output.
output = evalfis(fis,input,options) 
evaluates the fuzzy inference system using specified
evaluation options.
[output,fuzzifiedIn,ruleOut,aggregatedOut,ruleFiring] = evalfis( ___ ) 
returns intermediate results from the fuzzy inference process. 
This syntax is not supported when fis is a fistree object.

Evaluate Fuzzy Inference System
fis = readfis('tipper');
Evaluate the FIS when the first input is 2 and the second input is 1.
output = evalfis(fis,[2 1])
output = 7.0169

Evaluate FIS for Multiple Input Combinations
fis = readfis('tipper');
Specify the input combinations to evaluate using an array with one row per input combination.
input = [2 1;
         4 5;
         7 8];
output = evalfis(fis,input)
output = 3×1
7.0169
14.4585
20.3414
Each row of output is the defuzzified output value for the corresponding row of input.

Specify Number of Output Samples for FIS Evaluation
fis = readfis('tipper');
Create an evalfisOptions option set, specifying the number of samples in the output fuzzy sets.
options = evalfisOptions('NumSamplePoints',50);
output = evalfis(fis,[2 1],options);

Evaluate Tree of Fuzzy Inference Systems
fis1 = mamfis('Name','fis1','NumInputs',2,'NumOutputs',1);
fis2 = mamfis('Name','fis2','NumInputs',2,'NumOutputs',1);
con = ["fis1/output1" "fis2/input1"];
tree = fistree([fis1 fis2],con);

Create an evalfisOptions option set, specifying the number of samples in the output fuzzy sets.
options = evalfisOptions('NumSamplePoints',50);
y = evalfis(tree,[0.5 0.2 0.7],options)
y = 0.1553

Obtain Intermediate Fuzzy Inference Results
fis = readfis('tipper');
Evaluate the FIS, and return the intermediate inference results.
[output,fuzzifiedIn,ruleOut,aggregatedOut,ruleFiring] = evalfis(fis,[2 1]);
You can examine the intermediate results to understand or visualize the fuzzy inference process. 
For example, view the aggregated output fuzzy set, which is the fuzzy set that evalfis defuzzifies to find
the output value. 
Also, plot the defuzzified output value.
outputRange = linspace(fis.output.range(1),fis.output.range(2),length(aggregatedOut))';
plot(outputRange,aggregatedOut,[output output],[0 1])
xlabel('Tip')
ylabel('Output Membership')
legend('Aggregated output fuzzy set','Defuzzified output')
The length of aggregatedOutput corresponds to the number of sample points used to discretize
output fuzzy sets.


genfis
Generate fuzzy inference system object from data.

fis = genfis(inputData,outputData)
fis = genfis(inputData,outputData,options)

fis = genfis(inputData,outputData) 
returns a single-output Sugeno fuzzy inference system
(FIS) using a grid partition of the given input and output data.
fis = genfis(inputData,outputData,options) 
returns a FIS generated using the specified
input/output data and the options specified in options. 
You can generate fuzzy systems using grid
partitioning, subtractive clustering, or fuzzy c-means (FCM) clustering.

Generate Fuzzy Inference System Using Default Options
inputData = [rand(10,1) 10*rand(10,1)-5];
outputData = rand(10,1);
Generate a fuzzy inference system.
fis = genfis(inputData,outputData);
The generated system, fis, is created using grid partitioning with default options.

Generate FIS Using Grid Partitioning
inputData = [rand(10,1) 10*rand(10,1)-5];
outputData = rand(10,1);
Create a default genfisOptions option set for grid partitioning.
opt = genfisOptions('GridPartition');
Specify the following input membership functions for the generated FIS:
3 Gaussian membership functions for the first input variable
5 triangular membership functions for the second input variable
opt.NumMembershipFunctions = [3 5];
opt.InputMembershipFunctionType = ["gaussmf" "trimf"];
fis = genfis(inputData,outputData,opt);

Plot the input membership functions. 
Each input variable has the specified number and type of input
membership functions, evenly distributed over their input range.
[x,mf] = plotmf(fis,'input',1);
subplot(2,1,1)
plot(x,mf)
xlabel('input 1 (gaussmf)')
[x,mf] = plotmf(fis,'input',2);
subplot(2,1,2)
plot(x,mf)
xlabel('input 2 (trimf)')

Generate FIS Using Subtractive Clustering
load clusterDemo.dat
inputData = clusterDemo(:,1:2);
outputData = clusterDemo(:,3);
Create a genfisOptions option set and specify the range of influence for each data dimension.
Specify 0.5 and 0.25 as the range of influence for the first and second input variables. Specify 0.3
as the range of influence for the output data.
opt = genfisOptions('SubtractiveClustering',...
'ClusterInfluenceRange',[0.5 0.25 0.3]);
fis = genfis(inputData,outputData,opt);

The generated FIS contains one rule for each cluster.
showrule(fis)
ans = 3x83 char array
'1. If (in1 is in1cluster1) and (in2 is in2cluster1) then (out1 is out1cluster1) (1)'
'2. If (in1 is in1cluster2) and (in2 is in2cluster2) then (out1 is out1cluster2) (1)'
'3. If (in1 is in1cluster3) and (in2 is in2cluster3) then (out1 is out1cluster3) (1)'

Generate FIS Using FCM Clustering
load clusterDemo.dat
inputData = clusterDemo(:,1:2);
outputData = clusterDemo(:,3);
Create a genfisOptions option set for FCM Clustering, specifying a Mamdani FIS type.
opt = genfisOptions('FCMClustering','FISType','mamdani');
Specify the number of clusters.
opt.NumClusters = 3;
Suppress the display of iteration information to the Command Window.
opt.Verbose = 0;
fis = genfis(inputData,outputData,opt);

The generated FIS contains one rule for each cluster.
showrule(fis)
ans = 3x83 char array
'1. If (in1 is in1cluster1) and (in2 is in2cluster1) then (out1 is out1cluster1) (1)'
'2. If (in1 is in1cluster2) and (in2 is in2cluster2) then (out1 is out1cluster2) (1)'
'3. If (in1 is in1cluster3) and (in2 is in2cluster3) then (out1 is out1cluster3) (1)'

Plot the input and output membership functions.
[x,mf] = plotmf(fis,'input',1);
subplot(3,1,1)
plot(x,mf)
xlabel('Membership Functions for Input 1')
[x,mf] = plotmf(fis,'input',2);
subplot(3,1,2)
plot(x,mf)
xlabel('Membership Functions for Input 2')
[x,mf] = plotmf(fis,'output',1);
subplot(3,1,3)
plot(x,mf)
xlabel('Membership Functions for Output')


tunefis
Tune fuzzy inference system or tree of fuzzy inference systems.

fisOut = tunefis(fisIn,paramset,in,out)
fisOut = tunefis(fisIn,paramset,custcostfcn)
fisOut = tunefis( ___ ,options)
[fisOut,summary] = tunefis( ___ )

fisOut = tunefis(fisIn,paramset,in,out) 
tunes the fuzzy inference system fisin using the tunable parameter settings 
specified in paramset and the training data specified by in and out.
fisOut = tunefis(fisIn,paramset,custcostfcn) 
tunes the fuzzy inference system using a function handle to a custom cost function, custcostfcn.
fisOut = tunefis( ___ ,options) 
tunes the fuzzy inference system with additional options from
the object options created using tunefisOptions.
[fisOut,summary] = tunefis( ___ ) tunes the fuzzy inference system and returns additional
information about the tuning algorithm in summary.

Tune a Fuzzy Inference System.
x = (0:0.1:10)';
y = sin(2*x)./exp(x/5);
options = genfisOptions('GridPartition');
options.NumMembershipFunctions = 5;
fisin = genfis(x,y,options);

Obtain the tunable settings of inputs, outputs, and rules of the fuzzy inference system.
[in,out,rule] = getTunableSettings(fisin);

fisout = tunefis(fisin,[in;out],x,y,tunefisOptions("Method","anfis"));

Tune Specific Parameter Setting of Fuzzy Inference System.
x = (0:0.1:10)';
y = sin(2*x)./exp(x/5);
options = genfisOptions('GridPartition');
options.NumMembershipFunctions = 5;
fisin = genfis(x,y,options);
Obtain the tunable settings of inputs, outputs, and rules of the fuzzy inference system.
[in,out,rule] = getTunableSettings(fisin);
Tune the rule parameter only. In this example, the pattern search method is used.
fisout = tunefis(fisin,rule,x,y,tunefisOptions("Method","patternsearch"));

Learn Rules for FIS.
You can configure tunefis to learn the rules of a fuzzy system. 
For this example, learn rules for a tipping FIS.
fisin = readfis('tipper');
Generate training data using this FIS.
x = 10*rand(100,2);
y = evalfis(fisin,x);
Remove the rules from the FIS.
fisin.Rules = [];
To learn rules, set the OptimizationType option of tunefisOptions to "learning".
options = tunefisOptions( ...
"OptimizationType","learning", ...
"Display","none");
Set the maximum number of rules in the tuned FIS to 5.
options.NumMaxRules = 5;
Learn the rules without tuning any membership function parameters.
fisout = tunefis(fisin,[],x,y,options);

Tune a Fuzzy Inference System with Custom Parameter Settings.
x = (0:0.1:10)';
y = sin(2*x)./exp(x/5);
options = genfisOptions('GridPartition');
options.NumMembershipFunctions = 5;
fisin = genfis(x,y,options);
Obtain the tunable settings of inputs, outputs, and rules of the fuzzy inference system.
[in,out,rule] = getTunableSettings(fisin);
You can tune with custom parameter settings using setTunable or dot notation.

Do not tune input 1.
in(1) = setTunable(in(1),false);

For output 1:
do not tune membership functions 1 and 2,
do not tune membership function 3,
set the minimum parameter range of membership function 4 to -2,
set the maximum parameter range of membership function 5 to 2.
out(1).MembershipFunctions(1:2) = setTunable(out(1).MembershipFunctions(1:2),false);
out(1).MembershipFunctions(3).Parameters.Free = false;
out(1).MembershipFunctions(4).Parameters.Minimum = -2;
out(1).MembershipFunctions(5).Parameters.Maximum = 2;

For the rule settings,
do not tune rules 1 and 2,
set the antecedent of rule 3 to non-tunable,
allow NOT logic in the antecedent of rule 4,
do not ignore any outputs in rule 3.
rule(1:2) = setTunable(rule(1:2),false);
rule(3).Antecedent.Free = false;
rule(4).Antecedent.AllowNot = true;
rule(3).Consequent.AllowEmpty = false;

Set the maximum number of iterations to 20 and tune the fuzzy inference system.
opt = tunefisOptions("Method","particleswarm");
opt.MethodOptions.MaxIterations = 20;
fisout = tunefis(fisin,[in;out;rule],x,y,opt);

Prevent Overfitting Using K-Fold Cross-Validation.
To prevent the overfitting of your tuned FIS to your training data using k-fold cross validation.
This training data set has one input and one output.
load fuzex1trnData.dat
opt = genfisOptions('GridPartition');
opt.NumMembershipFunctions = 4;
opt.InputMembershipFunctionType = "gaussmf";
inputData = fuzex1trnData(:,1);
outputData = fuzex1trnData(:,2);
fis = genfis(inputData,outputData,opt);

For reproducibility, set the random number generator seed.
rng('default')
Configure the options for tuning the FIS. Use the default tuning method with a maximum of 30
iterations.
tuningOpt = tunefisOptions;
tuningOpt.MethodOptions.MaxGenerations = 30;

Configure the following options for using k-fold cross validation.
Use a k-fold value of 3.
Compute the moving average of the validation cost using a window of length 2.
Stop each training-validation iteration when the average cost is 5%
% greater than the current minimum cost.
tuningOpt.KFoldValue = 3;
tuningOpt.ValidationWindowSize = 2;
tuningOpt.ValidationTolerance = 0.05;

Obtain the settings for tuning the membership function parameters of the FIS.
[in,out] = getTunableSettings(fis);
[outputFIS,info] = tunefis(fis,[in;out],inputData,outputData,tuningOpt);

Evaluate the FIS for each of the training input values.
outputTuned = evalfis(outputFIS,inputData);

Plot the output of the tuned FIS along with the expected training output.
plot([outputData,outputTuned])
legend("Expected Output","Tuned Output","Location","southeast")
xlabel("Data Index")
ylabel("Output value")

Tune FIS Tree.
Create a FIS tree to model (sin x + cos x)/exp x. 

Create fis1 with two inputs, both with range [0, 10] and three MFs each. Use a smooth,
differentiable MF, such as gaussmf, to match the characteristics of the data type you are modeling.
fis1 = sugfis("Name","fis1");
fis1 = addInput(fis1,[0 10], ...
"NumMFs",3, ...
"MFType","gaussmf");
fis1 = addInput(fis1,[0 10], ...
"NumMFs",3, ...
"MFType","gaussmf");
Add an output with the range [–1.5, 1.5] having nine MFs corresponding to the nine possible input
MF combinations. Set the output range according to the possible values of sin x + cos x .
fis1 = addOutput(fis1,[-1.5 1.5],"NumMFs",9);

Create fis2 with two inputs. Set the range of the first input to [–1.5, 1.5], which matches the range
of the output of fis1. The second input is the same as the inputs of fis1. Therefore, use the same
input range, [0, 10]. Add three MFs for each of the inputs.
fis2 = sugfis("Name","fis2");
fis2 = addInput(fis2,[-1.5 1.5], ...
"NumMFs",3, ...
"MFType","gaussmf");
fis2 = addInput(fis2,[0 10], ...
"NumMFs",3, ...
"MFType","gaussmf");
Add an output with range [0, 1] and nine MFs. The output range is set according to the possible
values of (sin x + cos x)/exp x.
fis2 = addOutput(fis2,[0 1],"NumMFs",9);

Connect the inputs and the outputs as shown in the diagram. The first output of fis1 connects to
the first input of fis2. The inputs of fis1 connect to each other and the second input of fis1
connects to the second input of fis2.
con1 = ["fis1/output1" "fis2/input1"];
con2 = ["fis1/input1" "fis1/input2"];
con3 = ["fis1/input2" "fis2/input2"];

Create a FIS tree using the specified FISs and connections.
fisT = fistree([fis1 fis2],[con1;con2;con3]);
Add an additional output to the FIS tree to access the output of fis1.
fisT.Outputs = ["fis1/output1";fisT.Outputs];

For this example, generate input and output training data using the known mathematical operations.
Generate data for both the intermediate and final output of the FIS tree.
x = (0:0.1:10)';
y1 = sin(x)+cos(x);
y2 = y1./exp(x);
y = [y1 y2];

Learn the rules of the FIS tree using particle swarm optimization, which is a global optimization method.
options = tunefisOptions( ...
"Method","particleswarm", ...
"OptimizationType","learning");
This tuning step uses a small number of iterations to learn a rule base without overfitting the training data.
options.MethodOptions.MaxIterations = 5;
rng("default") % for reproducibility
fisTout1 = tunefis(fisT,[],x,y,options);

Tune all the FIS tree parameters at once using pattern search, which is a local optimization method.
options.Method = "patternsearch";
options.MethodOptions.MaxIterations = 25;

Use getTunableSettings to obtain input, output, and rule parameter settings from the FIS tree.
[in,out,rule] = getTunableSettings(fisTout1);
fisTout2 = tunefis(fisTout1,[in;out;rule],x,y,options);

The optimization cost is lower after the second tuning process.
Evaluate the FIS tree using the input training data.
yOut = evalfis(fisTout2,x);
Plot the final output along with the corresponding output training data.
plot(x,y(:,2),"-",x,yOut(:,2),"-")
legend("Training Data","FIS Tree Output")

The results do not perform well at the beginning and end of the input range. 
To improve performance, you could try:
Increasing the number of training iterations in each stage of the tuning process.
Increasing the number of membership functions for the input and output variables.
Using a custom cost function to model the known mathematical operations. 


Input Arguments

fisIn — Fuzzy inference system
mamfis object | sugfis object | mamfistype2 object | sugfistype2 object | fistree object
Fuzzy inference system, specified as one of the following objects.
mamfis object — Mamdani fuzzy inference system
sugfis object — Sugeno fuzzy inference system
mamfistype2 object — Type-2 Mamdani fuzzy inference system
sugfistype2 object — Type-2 Sugeno fuzzy inference system
fistree object — Tree of interconnected fuzzy inference systems

paramset — Tunable parameter settings
array
Tunable parameter settings, specified as an array of input, output, and rule parameter settings in the
input FIS. 
To obtain these parameter settings, use the getTunableSettings function with the input fisin.
paramset can be the input, output, or rule parameter settings, or any combination of these settings.

in — Input training data
matrix
Input training data, specified as an m-by-n matrix, where m is the total number of input datasets and
n is the number of inputs. The number of input and output datasets must be the same.

out — Output training data
matrix
Output training data, specified as an m-by-q matrix, where m is the total number of output datasets
and q is the number of outputs. The number of input and output datasets must be the same.

options — FIS tuning options
tunefisOptions object
FIS tuning options, specified as a tunefisOptions object. You can specify the tuning algorithm
method and other options for the tuning process.

custcostfcn — custom cost functions
function handle
Custom cost function, specified as a function handle. The custom cost function evaluates fisout to
calculate its cost with respect to an evaluation criterion, such as input/output data. 
custcostfcn must accept at least one input argument for fisout and returns a cost value. 
You can provide an anonymous function handle to attach additional data for cost calculation, as described in this example:
function fitness = custcost(fis,trainingData)
...
end
custcostfcn = @(fis)custcost(fis,trainingData);


Output Arguments

fisOut — Tuned fuzzy inference system
mamfis object | sugfis object | mamfistype2 object | sugfistype2 object | fistree object
Tuned fuzzy inference system, returned as one of the following objects.
mamfis object — Mamdani fuzzy inference system
sugfis object — Sugeno fuzzy inference system
mamfistype2 object — Type-2 Mamdani fuzzy inference system
sugfistype2 object — Type-2 Sugeno fuzzy inference system
fistree object — Tree of interconnected fuzzy inference systems
fisout is the same type of FIS as fisin.

summary — Tuning algorithm summary
structure
Tuning algorithm summary, specified as a structure containing the following fields:
tuningOutputs — Algorithm-specific tuning information
totalFunctionCount — Total number of evaluations of the optimization cost function
totalRuntime — Total execution time of the tuning process in seconds
errorMessage — Any error message generated when updating fisin with new parameter values
tuningOutputs is a structure that contains tuning information for the algorithm specified in options. 

The fields in tuningOutputs depend on the specified tuning algorithm.
When using k-fold cross validation:
tuningOutputs is an array of k structures, each containing the tuning information for one training-validation iteration.
totalFunctionCount and totalRuntime include the total function cost function evaluations
and total run time across all k training-validation iterations.





