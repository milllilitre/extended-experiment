% optimize PSO parameters using PSO

close all; clearvars; clc;
load('referenceMatrix.mat');
load('database.mat');
%% init
% constants
tagPopulation = 34;
rows = 1;
columns = 5;
channels = [920.625];
antennaNum = 2;
valNum = size(database, 3);
drawPics = 0;

tagLoc = zeros(tagPopulation,3);
tagLoc(1:17,1) = 0.6:0.15:3;
tagLoc(1:17,2) = 3.6 * ones(17,1);
tagLoc(1:17,3) = ones(17,1) * 0.43;
tagLoc(18:34,1) = 0 * ones(17,1);
tagLoc(18:34,2) = 0.6:0.15:3;
tagLoc(18:34,3) = 0.43 * ones(17,1);


readerLoc = [1.8 0 0.59
    3.6 1.8 0.59];
load('realLoc.mat');
% calculate negative phase
referenceMatrix(:,2) = -referenceMatrix(:,2);
database(:,2,:) = -database(:,2,:);

% remove nans
nanValue = zeros(2, valNum);
nanValue(:,:) = nanmin(database);
nanValue = min(transpose(nanValue));
nanLoc = isnan(database);
database(nanLoc) = 0;
nanLoc(:,2,:) = 0;
database(nanLoc) = nanValue(1);

% transform reference matrix into vector form
refVecMatrix = zeros(tagPopulation, 2); % reference Matrix, in signal vector form
[refVecMatrix(:, 1), refVecMatrix(:, 2)] = rssiPhase2IQ(referenceMatrix(:, 1), referenceMatrix(:, 2));

%% process data
% first, all data divide by 2 (which means sqrt for channel parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% referenceMatrix = referenceMatrix / 2;
% database = database / 2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for i = 5:5
%     figure(); subplot(2,1,1);plot(database(:,1,i), '*'); title('raw phase measurements');
%     subplot(2,1,2); plot(referenceMatrix(:,2), '*'); title('reference phase measurements');
% end


% then, all data minus reference signal (in vector form)
removeRefMatrix = zeros(size(database));
tmpMat = zeros(size(referenceMatrix));
removeRefMethod = 1;
if(removeRefMethod == 1)
    for fileIndex = 1:valNum
        tagVecMatrix = zeros(tagPopulation, 2); % reference Matrix, in signal vector form
        [tagVecMatrix(:, 1), tagVecMatrix(:, 2)] = rssiPhase2IQ(database(:, 1, fileIndex), database(:, 2, fileIndex));
        tagVecMatrix = tagVecMatrix - refVecMatrix;
        [removeRefMatrix(:,1,fileIndex), removeRefMatrix(:,2,fileIndex)] = ...
            IQ2rssiPhase(tagVecMatrix(:,1), tagVecMatrix(:,2));
    end
else
    % wrong!! It is reference matrix that need to add pi, not tag
    % measurement matrix!!!
    for fileIndex = 1:valNum
        removeRefMatrix(:,:,fileIndex) = database(:,:,fileIndex);
        removeRefMatrix(:,2,fileIndex) = removeRefMatrix(:,2,fileIndex) + pi;
        [removeRefMatrix(:,1,fileIndex), removeRefMatrix(:,2,fileIndex)] = ...
            vecSum(removeRefMatrix(:,1,fileIndex), removeRefMatrix(:,2,fileIndex), referenceMatrix(:,1), referenceMatrix(:,2));
    end
end


if(drawPics)
    for i = 1:valNum
    % for i = 5:5
        figure(); plot(removeRefMatrix(:,2,i), '*'); title(['remove reference ' num2str(i)]);
    end
end
% calculate difference of adjacent tag pairs,注意，1~17标签一组，18~34标签一组
tagDeltaMatrix = zeros(tagPopulation - 2, 2, valNum);
tagDeltaMatrix(1:16, :, :) = removeRefMatrix(2:17, :, :) - removeRefMatrix(1:16, :, :);
tagDeltaMatrix(17:32, :, :) = removeRefMatrix(19:34, :, :) - removeRefMatrix(18:33, :, :);

% for i = 5:5
%     figure(); plot(tagDeltaMatrix(:,2,i), '*'); title('dtag, without regulation');
% end



%% remove difference caused by tag to reader link(phase and amp)
if(1)
    Drt = zeros(tagPopulation, 1);
    for tagIndex = 1:tagPopulation
        if(tagIndex <= 17)
            Drt(tagIndex) = norm(tagLoc(tagIndex,:) - readerLoc(1,:));
        else
            Drt(tagIndex) = norm(tagLoc(tagIndex,:) - readerLoc(2,:));
        end
    end
    DeltaDrt = zeros(tagPopulation - 2, 1);
    DeltaDrt(1:16) = Drt(2:17) - Drt(1:16);
    DeltaDrt(17:32) = Drt(19:34) - Drt(18:33);
    for fileIndex = 1:valNum
        tagDeltaMatrix(:, 2, fileIndex) = mod(tagDeltaMatrix(:, 2, fileIndex) - DeltaDrt, 2 * pi);  % phase
    end
    % should not square root directly, this is rssi value, should transform
end

% for i = 5;5
%     figure(); plot(tagDeltaMatrix(:,2,i), '*'); title('remove tag to reader diff');
% end


%% phase resolve according to certain priori knowledge
% actual phase of a line of tag should be a 'V' like shape
% first, find the v-zone of the phase: use a sliding window, if the rounded
% distance between the 3 points are shortest, then select the center point
% as the point of V-shape
if(1)
    phaseMatrix = zeros(size(removeRefMatrix));
    windowSize = 4;
    windowMatrix = zeros(tagPopulation - 2 * windowSize, valNum);
    centerPosi = zeros(valNum, 2);
    
    % regulate tagDeltaMatrix
    tmpMatrix = tagDeltaMatrix(:,:,:) > pi;    
    tmpMatrix(:,1,:) = 0;
    tagDeltaMatrix(tmpMatrix) = tagDeltaMatrix(tmpMatrix) - 2 * pi;
    tmpMatrix = tagDeltaMatrix(:,2,:) < -pi;
    tmpMatrix(:,1,:) = 0;
    tagDeltaMatrix(tmpMatrix) = tagDeltaMatrix(tmpMatrix) + 2 * pi;
    
    % calculate phaseMatrix according to tagDeltaMatrix, assume the first
    % phase value is 0;
    for fileIndex = 1:valNum
        for tagIndex = 1:16
            phaseMatrix(tagIndex + 1, 2, fileIndex) = phaseMatrix(tagIndex, 2, fileIndex) + tagDeltaMatrix(tagIndex, 2, fileIndex);
        end
        for tagIndex = 17:32
            phaseMatrix(tagIndex + 2, 2, fileIndex) = phaseMatrix(tagIndex + 1, 2, fileIndex) + tagDeltaMatrix(tagIndex, 2, fileIndex);
        end
    end
    
    
    % find centerPosi
    for fileIndex = 1:valNum
        for tagIndex = 1:(17-windowSize)
            windowMatrix(tagIndex, fileIndex) = sum(abs(tagDeltaMatrix(tagIndex:(tagIndex + windowSize - 1), 2, fileIndex)));
        end
        [sorted, indices] = sort(windowMatrix(1:(17 - windowSize), fileIndex));
        centerPosi(fileIndex, 1) = indices(1) + floor(windowSize / 2);
        
        for tagIndex = (18 - windowSize):(tagPopulation - 2 * windowSize)
            windowMatrix(tagIndex, fileIndex) = sum(abs(tagDeltaMatrix((tagIndex + windowSize - 1):(tagIndex + 2 * windowSize - 2), 2, fileIndex)));
        end
        [sorted, indices] = sort(windowMatrix((18-windowSize):(tagPopulation - 2 * windowSize), fileIndex));
        centerPosi(fileIndex, 2) = indices(1) + floor(windowSize / 2) + 17 ;
    end
    
    % regulate phase according to centerPosi using threshold methods
    for fileIndex = 1:valNum
        tmpVec = transpose(phaseMatrix(centerPosi(fileIndex, 1):-1:1, 2, fileIndex));
        tmpVec = phaseRegulateVec(tmpVec);
        phaseMatrix(centerPosi(fileIndex, 1):-1:1, 2, fileIndex) = transpose(tmpVec);
        
        tmpVec = transpose(phaseMatrix(centerPosi(fileIndex, 1):17, 2, fileIndex));
        tmpVec = phaseRegulateVec(tmpVec);
        phaseMatrix(centerPosi(fileIndex, 1):17, 2, fileIndex) = transpose(tmpVec);
        
        tmpVec = transpose(phaseMatrix(centerPosi(fileIndex, 2):-1:18, 2, fileIndex));
        tmpVec = phaseRegulateVec(tmpVec);
        phaseMatrix(centerPosi(fileIndex, 2):-1:18, 2, fileIndex) = transpose(tmpVec);
        
        tmpVec = transpose(phaseMatrix(centerPosi(fileIndex, 2):34, 2, fileIndex));
        tmpVec = phaseRegulateVec(tmpVec);
        phaseMatrix(centerPosi(fileIndex, 2):34, 2, fileIndex) = transpose(tmpVec);
    end
    
    if(drawPics)
        for fileIndex = 1:valNum
            figure(); plot(phaseMatrix(:,2,fileIndex), '*'); title(['regulate phase ' num2str(fileIndex)]);
            hold on;
            plot(centerPosi(fileIndex, 1:2), phaseMatrix(centerPosi(fileIndex, 1:2), 2, fileIndex), 'o');
        end
    end
end

    % regulate phase according to centerPosi using modified pool adjacent
    % violators algorithm : that is, if two phase values 
if(1)
    for fileIndex = 1:valNum
        tmpVec = transpose(phaseMatrix(centerPosi(fileIndex, 1):-1:1, 2, fileIndex));
        tmpVec = lsqisotonic(1:centerPosi(fileIndex, 1), tmpVec);
        phaseMatrix(centerPosi(fileIndex, 1):-1:1, 2, fileIndex) = transpose(tmpVec);
        
        tmpVec = transpose(phaseMatrix(centerPosi(fileIndex, 1):17, 2, fileIndex));
        tmpVec = lsqisotonic(centerPosi(fileIndex, 1):17, tmpVec);
        phaseMatrix(centerPosi(fileIndex, 1):17, 2, fileIndex) = transpose(tmpVec);
        
        tmpVec = transpose(phaseMatrix(centerPosi(fileIndex, 2):-1:18, 2, fileIndex));
        tmpVec = lsqisotonic(18:centerPosi(fileIndex, 2), tmpVec);
        phaseMatrix(centerPosi(fileIndex, 2):-1:18, 2, fileIndex) = transpose(tmpVec);
        
        tmpVec = transpose(phaseMatrix(centerPosi(fileIndex, 2):34, 2, fileIndex));
        tmpVec = lsqisotonic(centerPosi(fileIndex, 2):34, tmpVec);
        phaseMatrix(centerPosi(fileIndex, 2):34, 2, fileIndex) = transpose(tmpVec);
    end
    
    if(drawPics)
        for i = 1:valNum
            figure(); plot(phaseMatrix(:,2,i), '*'); title(['after isotonic regression ' num2str(i)]);
        end
    end
end



%% optimize PSO定位 
if(1)
    methods = 1;
%     ------给定初始化条件----------------------------------------------  
%     learnRate1 = 1.4962 / 2;             %加速常数即学习因子1  
%     learnRate1 = 2;
%     learnRate2 = 2;             %加速常数即学习因子2  
    learnRate1 = 0.4962;
    learnRate2 = 0.4962;
    inertia = 0.3298;              %惯性权重  
%     inertia = 0.5;
    maxIter = 100;           %最大迭代次数  
    dimention = 3;                  %搜索空间维数（测试函数fitnessFcn中未知数个数）  
    population = 30;                  %初始化群体个体数目  
    accuracy = 0.2;           %设置精度(在已知最小值时候用)  
    upperBound = [3.6 3.6 1];
    lowerBound = [0 0 0];
    rawResult = zeros(valNum, 2, methods);
    errorMat = zeros(valNum, methods);
    avaDists = zeros(methods);
    timeVec = zeros(methods, 1);
    convergMaxIter = 100;
    iters = zeros(valNum, methods);
    
    % ---- commence PSO optim ----
    data = transpose(tagDeltaMatrix(:, 2, fileIndex));
    [Pbest, result] = psoPsoSimple2(data, learnRate1, ...
        learnRate2, inertia, maxIter, dimention, population, ...
        accuracy, upperBound, lowerBound, convergMaxIter);
    result
    % show cummu prab comparison and average localization error of each
    % method
end



%% tagDeltaMatrix tranform into amp multiplier and distance
if(0)
    %??????????????????????????????????????????????????????
    channelIndex = 1;
    waveLength = 3e8 / (channels(channelIndex) * 1e6);
    distMat = zeros(size(tagDeltaMatrix));
    % rssi -> amp multiplier
    distMat(:,1,:) = 10 .^ (tagDeltaMatrix(:,1,:) / 10);
    % phase -> distance
    distMat(:,2,:) = tagDeltaMatrix(:,2,:) / 2 / pi * waveLength;
    %??????????????????????????????????????????????????????
end




