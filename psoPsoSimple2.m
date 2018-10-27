% 创建并初始化一个n维的粒子群, optimized using matrix form
% 20180624 added a mechenism to weight the fitnessFcn with distance

function [Pbest, result, iters] = psoPsoSimple2(data, learnRate1, learnRate2, ...
    inertia, maxIter, dimention, population, accuracy, upperBound, lowerBound,...
    convergMaxIter)
load('tagDeltaMatrix.mat');
data = tagDeltaMatrix;
% repeat  
%      for 每个粒子i=1,2,…n do  
%            //设置个体最佳位置  
%            if f(i)<y then  
%                 y=f(i);  
%            end  
%            //设置全局最佳位置  
%            if y<Y then  
%                 Y=y;  
%            end  
%      end  
%      for 每个粒子i=1,2,…n do  
%          用速度方程更新速度  
%          用位置方程更新位置  
%      end  
% until 满足终止条件  
% gbest PSO的matlab实现代码如下：  
%------初始格式化--------------------------------------------------  
format long;  
%------给定初始化条件----------------------------------------------  
c1 = learnRate1;
c2 = learnRate2;
w = inertia;
MaxDT = maxIter;
D = dimention;
N = population;
eps = accuracy;

%------初始化种群的个体(可以在这里限定位置和速度的范围)------------  
x = randn(N, D) + 1.8; %产生一个服从正态分布的随机数作为初始化位置  
x(1,:) = [1.4962 1.4962 0.7298];

x = checkBounds(x, upperBound, lowerBound);
v = randn(N, D) + 1.8; %产生一个服从正态分布的随机数作为初始化速度  
v = v / 10;
%------先计算各个粒子的适应度，并初始化个体最优位置y和全局最优位置Pg--------  
%计算适应度，测试函数为fitnessFcn, calculate using matrix form
p = fitnessFcnMat(x, data); % p records the best fitness value of each particle
y = x;  %初始化个体最优位置y为在时间步t=0时的粒子位置  
[tmpVal, tmpIndex] = min(transpose(p));
Pg = y(tmpIndex,:);             %Pg为全局最优位置  
Pbest = zeros(1, maxIter);
Pbest(1) = tmpVal;
convergIter = 0;

%------进入主要循环，按照公式依次迭代，直到满足精度要求 or reach max iters------------  
for t = 1:MaxDT
    v = v * w + c1 * rand(N,D) .* (y - x) + c2 * rand(N, D) .* (repmat(Pg, N, 1) - x);
    x = x + v;
    [x, v] = checkBoundsBouncingWall(x, v, upperBound, lowerBound);
    %更新适应度 and best position of each particle
    tmpP = fitnessFcnMat(x, data);
    tmpComp = tmpP < p;
    p(tmpComp) = tmpP(tmpComp); % refresh best fitness of each particle
    tmpComp = repmat(tmpComp, 1, D);
    y(tmpComp) = x(tmpComp);    % refresh best position of each particle
    %更新群体最佳位置  
    [tmpVal, tmpIndex] = min(transpose(p));
    Pg = y(tmpIndex, :);        % find location of current best particle
    if(t > 1)   % if no group best is found, update convergIter flag
        if(tmpVal == Pbest(t - 1))
            convergIter = convergIter + 1;
        else
            convergIter = 0;
        end
    else
        if(tmpVal == Pbest(1))
            convergIter = convergIter + 1;
        else
            convergIter = 0;
        end
    end
    Pbest(t) = p(tmpIndex);
    % if group best fitness already reached accuracy, exit loop
    if(Pbest(t) < accuracy)
        break;
    end
    % if already some iters without updating pbest, exit loop
    if(convergIter >= convergMaxIter)
        break;
    end
end
iters = t;
Pbest = Pbest(1:t);
result = Pg;
end

function outputMat = fitnessFcnMat(posiMat, measurement)
load('realLoc.mat');
tagDeltaMatrix = measurement;
particleNum = size(posiMat, 1);
valNum = size(tagDeltaMatrix, 3);
learnRate1 = posiMat(1);
learnRate2 = posiMat(2);
inertia = posiMat(3);              %惯性权重  
%     inertia = 0.5;
maxIter = 50;           %最大迭代次数  
dimention = 2;                  %搜索空间维数（测试函数fitnessFcn中未知数个数）  
population = 30;                  %初始化群体个体数目  
accuracy = 1;           %设置精度(在已知最小值时候用)  
upperBound = [3.6 3.6];
lowerBound = [0 0];
rawResult = zeros(valNum, 2);
convergMaxIter = 5;
% iters = zeros(valNum, methods);
outputMat = zeros(particleNum, 1);
% ---- commence PSO optim ----
for particleIndex = 1:particleNum
    rawResult = zeros(valNum, dimention);
    for fileIndex = 1:valNum
        data = transpose(tagDeltaMatrix(:, 2, fileIndex));
        [Pbest, result] = psoSimple2(data, learnRate1, ...
            learnRate2, inertia, maxIter, dimention, population, ...
            accuracy, upperBound, lowerBound, convergMaxIter);
    %         figure(); plot(Pbest);
        rawResult(fileIndex, :) = result;
    end
    % calculate and show accuracy
    [distVec, avaDist, maxDist] = calculateAccuracy(rawResult(:,:), realLoc(:,1:2));
    outputMat(particleIndex) = avaDist;
end
end



function output = fitnessFcn(posi, measurement)
tagPopulation = 34;
tagLoc = zeros(tagPopulation,3);
tagLoc(1:17,1) = 0.6:0.15:3;
tagLoc(1:17,2) = 3.6 * ones(17,1);
tagLoc(1:17,3) = ones(17,1) * 0.43;
tagLoc(18:34,1) = 0 * ones(17,1);
tagLoc(18:34,2) = 0.6:0.15:3;
tagLoc(18:34,3) = 0.43 * ones(17,1);
readerLoc = [1.8 0 0.59
    3.6 1.8 0.59];
currentPosi = [posi 0.59];

tagPopulation = 34;
waveLength = 3e8 / 920.625e6;
deltaPhase = zeros(1, tagPopulation - 1);
for tagIndex = 1:16
    d1 = norm(tagLoc(tagIndex, :) - currentPosi);
    d2 = norm(tagLoc(tagIndex + 1, :) - currentPosi);
    deltaD = d2 - d1;
    deltaPhase(tagIndex) = deltaD / waveLength * 2 * pi;
end
for tagIndex = 17:(tagPopulation - 2)
    d1 = norm(tagLoc(tagIndex + 1, :) - currentPosi);
    d2 = norm(tagLoc(tagIndex + 2, :) - currentPosi);
    deltaD = d2 - d1;
    deltaPhase(tagIndex) = deltaD / waveLength * 2 * pi;
end
output = norm(measurement - deltaPhase);

end

