% input is rows of readings, and rows of loaction
function [resultLabel] = KNNlocFunc(k, trainData, trainLoc, valData)
    dimension = size(trainLoc, 2);    % dimention of localization
    nFeature = size(trainData, 2);
    trainDataNum = size(trainData, 1);
    valDataNum = size(valData, 1);

    %% main logic
    % training process: calculate average, here we asume all the training data
    % of a certain label are piled together
    %% normalization
    for i = 1:nFeature
        tmpTrainMax = max(trainData(:,i));
        tmpTrainMin = min(trainData(:,i));
        trainData(:,i) = (trainData(:,i) - tmpTrainMin) / (tmpTrainMax - tmpTrainMin);
        valData(:,i) = (valData(:,i) - tmpTrainMin) / (tmpTrainMax - tmpTrainMin);
    end
    
    
    if(0)
        i = 1;
        j = 1;
        currentLabel = trainLoc(1,:);
        currentSampleSum = zeros(1, nFeature);
        currentSampleNum = 0;
        dataSet = zeros(size(trainData));
        labelSet = zeros(size(trainLoc));
        labelSet(1,:) = trainLoc(1,:);
        dataSetLen = 1;
        while(i <= trainDataNum)
            if(trainLoc(i,:) == currentLabel)
                currentSampleSum = currentSampleSum + trainData(i, :);
                currentSampleNum = currentSampleNum + 1;
            else
                dataSet(j, :) = currentSampleSum / currentSampleNum;
                currentLabel = trainLoc(i,:);   % reset current label
                currentSampleNum = 0;
                j = j + 1;
                dataSetLen = dataSetLen + 1;
                labelSet(j, :) = trainLoc(i, :);
                currentSampleSum = trainData(i,:);
            end
            i = i + 1;
        end
        if(i == trainDataNum)
            dataSet(j,:) = currentSampleSum / currentSampleNum;
        end
        dataSet = dataSet(1:dataSetLen, :);
    %     size(dataSet)
        labelSet = labelSet(1:dataSetLen, :);
    %     size(labelSet)
    end
    resultLabel = zeros(valDataNum, dimension);
    for i = 1:valDataNum
%          resultLabel(i,:) = KNNlocSingle(k, dataSet, labelSet, valData(i,:));
       resultLabel(i,:) = KNNlocSingle(k, trainData, trainLoc, valData(i,:));
    end
end

function [resultLabel] = KNNlocSingle(k, dataSet, labelSet, valData)
    %% 
    %   inx 为 输入测试数据，data为样本数据，labels为样本标签
    data = dataSet;
    labels = labelSet;
    inx = valData;
    %%

    [datarow , datacol] = size(data);
    diffMat = repmat(inx,[datarow,1]);
    diffMat = diffMat - data ;
    distanceMat = sqrt(sum(diffMat.^2,2));
    [B , IX] = sort(distanceMat,'ascend');
%    len = min(k,length(B));
    len = k;
    resultLabel = vote(labels(IX(1:len),:));
%      resultLabel = cell2mat(mode(num2cell(resultLabel, 1)));
%     resultLabel = mode(resultLabel);
%     resultLabel = labels(IX(1:len),:)
end