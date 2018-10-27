close all; clearvars; clc;
%这个文件用来处理双墙面变高度50标签反射定位试验的数据，use PSO method
% !!!! phase value calculated using impinj R420 reader lowers when distance
% increases, which is the opposite of what it should be!!!
% according to formula, removed Hxt1/Hxt2 component
load('referenceMatrixVec.mat');
load('database.mat');   % tagPopulation * Rssi & phase * valNum * testNum
load('realLoc.mat');
%% init
% constants
tagPopulation = 50;
tagEachArray = 25;
% rows = 1;
% columns = 5;
channels = [920.625];
antennaNum = 2;
valNum = size(databaseFull, 3); 
testNum = 3;
antennas = 2;
lines_used_in_tiles = 120;
% toggle switches
drawPics = 0;

tagLocAll = zeros(tagPopulation, 3, testNum);
readerLocAll = zeros(antennaNum, 3, testNum);
tagLocAll(1:tagEachArray, 1, :) = repmat(-0.3:0.15:3.3, 1, 1, 3);
tagLocAll(1:tagEachArray, 2, :) = repmat(3.6 , tagEachArray, 1, 3);
tagLocAll((tagEachArray + 1):tagPopulation, 1, :) = repmat(3.6, tagEachArray, 1, 3);
tagLocAll((tagEachArray + 1):tagPopulation, 2, :) = repmat(transpose(-0.3:0.15:3.3), 1, 1, 3); 
tagLocAll(1:50, 3, 1) = ones(tagPopulation, 1, 1) * 1.3;
tagLocAll(1:50, 3, 2) = ones(tagPopulation, 1, 1) * 1.6;
tagLocAll(1:50, 3, 3) = ones(tagPopulation, 1, 1) * 1.85;  % 注意标签和阅读器高度不一

tagLocAll(28, 3, 1) = 1.1;
tagLocAll(27, 3, 2) = 1.6;

readerLocAll(:, :, 1) = [1.8 -0.1 1.3
    0 1.8 1.3];
readerLocAll(:, :, 2) = [1.8 -0.1 1.6
    0 1.8 1.6];
readerLocAll(:, :, 3) = [1.8 -0.1 1.9
    0 1.8 1.9];

% 阅读器位置[1.8 -0.1 1.3    0 1.8 1.3];

% calculate negative phase. -- because of the property of Impinj readers
referenceMatrixVec(:, 2, :) = - referenceMatrixVec(:, 2, :);
databaseFull(:, 2, :) = - databaseFull(:, 2, :);

% initialize some of the parameters and variables
methods = 4;
x_error_vec1 = zeros(1, testNum);
y_error_vec1 = zeros(1, testNum);
avaDists = zeros(methods, testNum);
medDists = zeros(methods, testNum);
variance = zeros(methods, testNum);
timeMatrix = zeros(methods, testNum);
x_error_vec2 = zeros(1, testNum);
y_error_vec2 = zeros(1, testNum);


for test_index = 1:testNum
    % initiate test-specific data
    tagLoc = zeros(tagPopulation, 3);
    tagLoc(:,:) = tagLocAll(:, :, test_index);
    readerLoc = zeros(2, 3);
    readerLoc(:, :) = readerLocAll(:, :, test_index);
    realLoc = zeros(valNum, 3);
    realLoc(:, :) = realLocAll(:, :, test_index);
    database = zeros(tagPopulation, 2, valNum);
    database(:, :, :) = databaseFull(:, :, :, test_index);
    referenceMatrix = zeros(tagPopulation, 2);
    referenceMatrix(:, :) = referenceMatrixVec(:, :, test_index + 2);  % 由于第一组比较特殊，测试了三个ref mat
    
    % remove nans
    nanValue = zeros(2, valNum);
    nanValue(:, :) = nanmin(database);
    nanValue = min(transpose(nanValue));
    nanLoc = isnan(database);
    database(nanLoc) = 0;
    nanLoc(:, 2, :) = 0;
    database(nanLoc) = nanValue(1);

    % transform reference matrix into vector form
    refVecMatrix = zeros(tagPopulation, 2); % reference Matrix, in signal vector form
    [refVecMatrix(:, 1), refVecMatrix(:, 2)] = rssiPhase2IQ(referenceMatrix(:, 1), referenceMatrix(:, 2));

    %% process data


    % all data minus reference signal (in vector form)
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

    % the first set of data is special because it has three ref mats, each
    % covers a protion of the data. So removeRefMatrix generation is a bit
    % more complex
    if (test_index == 1)
        
    end

    if(drawPics)
        for i = 1:valNum
        % for i = 5:5
            figure(); plot(removeRefMatrix(:,2,i), '*'); title(['remove reference ' num2str(i)]);
        end
    end
    % calculate difference of adjacent tag pairs,注意，1~17标签一组，18~34标签一组
    tagDeltaMatrix = zeros(tagPopulation - 2, 2, valNum);
    tagDeltaMatrix(1:(tagEachArray - 1), :, :) = removeRefMatrix(2:tagEachArray, :, :) - removeRefMatrix(1:(tagEachArray - 1), :, :);
    tagDeltaMatrix(tagEachArray:(tagPopulation - 2), :, :) = removeRefMatrix((tagEachArray + 2):tagPopulation, :, :) - removeRefMatrix((tagEachArray + 1):(tagPopulation - 1), :, :);


    %% remove difference caused by tag to reader link(phase and amp)
    if(1)
        Drt = zeros(tagPopulation, 1);
        for tagIndex = 1:tagPopulation
            if(tagIndex <= tagEachArray)
                Drt(tagIndex) = norm(tagLoc(tagIndex,:) - readerLoc(1,:));
            else
                Drt(tagIndex) = norm(tagLoc(tagIndex,:) - readerLoc(2,:));
            end
        end
        DeltaDrt = zeros(tagPopulation - 2, 1);
        DeltaDrt(1:(tagEachArray - 1)) = Drt(2:tagEachArray) - Drt(1:(tagEachArray - 1));
        DeltaDrt(tagEachArray:(tagPopulation - 2)) = Drt((tagEachArray + 2):tagPopulation) - Drt((tagEachArray + 1):(tagPopulation - 1));
        for fileIndex = 1:valNum
            tagDeltaMatrix(:, 2, fileIndex) = mod(tagDeltaMatrix(:, 2, fileIndex) - DeltaDrt, 2 * pi);  % phase
            tagDeltaMatrix(1:(tagEachArray - 1), 1, fileIndex) = tagDeltaMatrix(1:(tagEachArray - 1), 1, fileIndex) ...
                ./ Drt(2:tagEachArray) .* Drt(1:(tagEachArray - 1));  % amp
            tagDeltaMatrix(tagEachArray:(tagPopulation - 2), 1, fileIndex) = tagDeltaMatrix(tagEachArray:(tagPopulation - 2), 1, fileIndex) ...
                ./ Drt((tagEachArray + 2):tagPopulation) .* Drt((tagEachArray + 1):(tagPopulation - 1));
        end
    end


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
            for tagIndex = 1:(tagEachArray - 1)
                phaseMatrix(tagIndex + 1, 2, fileIndex) = ...
                    phaseMatrix(tagIndex, 2, fileIndex) + tagDeltaMatrix(tagIndex, 2, fileIndex);
            end
            for tagIndex = tagEachArray:(tagPopulation - 2)
                phaseMatrix(tagIndex + 2, 2, fileIndex) = ...
                    phaseMatrix(tagIndex + 1, 2, fileIndex) + tagDeltaMatrix(tagIndex, 2, fileIndex);
            end
        end


        % find centerPosi
        for fileIndex = 1:valNum
            for tagIndex = 1:(tagEachArray - windowSize)
                windowMatrix(tagIndex, fileIndex) = ...
                    sum(abs(tagDeltaMatrix(tagIndex:(tagIndex + windowSize - 1), 2, fileIndex)));
            end
            [sorted, indices] = sort(windowMatrix(1:(tagEachArray - windowSize), fileIndex));
            centerPosi(fileIndex, 1) = indices(1) + floor(windowSize / 2);

            for tagIndex = (tagEachArray + 1 - windowSize):(tagPopulation - 2 * windowSize)
                windowMatrix(tagIndex, fileIndex) = ...
                    sum(abs(tagDeltaMatrix((tagIndex + windowSize - 1):(tagIndex + 2 * windowSize - 2), 2, fileIndex)));
            end
            [sorted, indices] = sort(windowMatrix((tagEachArray + 1 - windowSize):(tagPopulation - 2 * windowSize), fileIndex));
            centerPosi(fileIndex, 2) = indices(1) + floor(windowSize / 2) + tagEachArray ;
        end

        % regulate phase according to centerPosi using threshold methods
        for fileIndex = 1:valNum
            tmpVec = transpose(phaseMatrix(centerPosi(fileIndex, 1):-1:1, 2, fileIndex));
            tmpVec = phaseRegulateVec(tmpVec);
            phaseMatrix(centerPosi(fileIndex, 1):-1:1, 2, fileIndex) = transpose(tmpVec);

            tmpVec = transpose(phaseMatrix(centerPosi(fileIndex, 1):tagEachArray, 2, fileIndex));
            tmpVec = phaseRegulateVec(tmpVec);
            phaseMatrix(centerPosi(fileIndex, 1):tagEachArray, 2, fileIndex) = transpose(tmpVec);

            tmpVec = transpose(phaseMatrix(centerPosi(fileIndex, 2):-1:(tagEachArray + 1), 2, fileIndex));
            tmpVec = phaseRegulateVec(tmpVec);
            phaseMatrix(centerPosi(fileIndex, 2):-1:(tagEachArray + 1), 2, fileIndex) = transpose(tmpVec);

            tmpVec = transpose(phaseMatrix(centerPosi(fileIndex, 2):tagPopulation, 2, fileIndex));
            tmpVec = phaseRegulateVec(tmpVec);
            phaseMatrix(centerPosi(fileIndex, 2):tagPopulation, 2, fileIndex) = transpose(tmpVec);
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

            tmpVec = transpose(phaseMatrix(centerPosi(fileIndex, 1):tagEachArray, 2, fileIndex));
            tmpVec = lsqisotonic(centerPosi(fileIndex, 1):tagEachArray, tmpVec);
            phaseMatrix(centerPosi(fileIndex, 1):tagEachArray, 2, fileIndex) = transpose(tmpVec);

            tmpVec = transpose(phaseMatrix(centerPosi(fileIndex, 2):-1:(tagEachArray + 1), 2, fileIndex));
            tmpVec = lsqisotonic((tagEachArray + 1):centerPosi(fileIndex, 2), tmpVec);
            phaseMatrix(centerPosi(fileIndex, 2):-1:(tagEachArray + 1), 2, fileIndex) = transpose(tmpVec);

            tmpVec = transpose(phaseMatrix(centerPosi(fileIndex, 2):tagPopulation, 2, fileIndex));
            tmpVec = lsqisotonic(centerPosi(fileIndex, 2):tagPopulation, tmpVec);
            phaseMatrix(centerPosi(fileIndex, 2):tagPopulation, 2, fileIndex) = transpose(tmpVec);
        end

        if(drawPics)
            for i = 1:valNum
                figure(); plot(phaseMatrix(:,2,i), '*'); title(['after isotonic regression ' num2str(i)]);
            end
        end
    end



    %% PSO定位
    if(1)
    %     ------给定初始化条件----------------------------------------------  
    if(1)   %PSO optimized parameter
        learnRate1 = 0.918746094077224;
        learnRate2 = 1.305432218802911;
        inertia = 0.878;
    end
    if(0)   % according to <<Handbook of SI: Concepts, Principles and Applications>>
        % but result is worse than PSO optimized parameters
        learnRate1 = 0.5 + log(2);
        learnRate2 = learnRate1;
        inertia = 1 / 2 / log(2);
    end

    %     inertia = 0.5;
        maxIter = 20;           %最大迭代次数 
        dimention = 2;                  %搜索空间维数（测试函数fitnessFcn中未知数个数）  
        population = 13;                  %初始化群体个体数目  
        accuracy = 1;           %设置精度(在已知最小值时候用)  
        upperBound = [3.6 3.6];
        lowerBound = [0 0];
        rawResult = zeros(valNum, 2, methods);
        errorMat = zeros(valNum, methods);
        timeVec = zeros(methods, 1);
        convergMaxIter = 20;
        iters = zeros(valNum, methods);
        % --------------------initialize LSQ optim method------------------
        fun1 = @method2Func1;                   % fitness used to locate x
        fun2 = @method2Func2;                   % fitness used to locate y
        fun3 = @method2Func;                    % fitness used to locate both x and y
%         xdata = zeros(6,(tagPopulation - 2));   % all the known parameters
        xdata = tagLoc;
        xdata = [xdata
            realLoc(1,:)];
        ydata = zeros((tagPopulation - 2), 1);  % c
        ub = [3.6, 3.6, pi * 2];
        lb = [0, 0, 0];
    %     x0 = transpose([1.8 1.8 0]);
    %     x0 = transpose([0, 0, 0]);
        options = optimoptions('lsqcurvefit',...
            'Algorithm','levenberg-marquardt',...
            'MaxIterations', 20,...
            'OptimalityTolerance', 1e-3,...
            'FunctionTolerance',1e-3, ...
            'Display', 'off');

        % ---- commence localization using different methods ----
        result_x1 = zeros(valNum, 2);
        result_y1 = zeros(valNum, 2);
        result_x2 = zeros(valNum, 2);
        result_y2 = zeros(valNum, 2);
        for methodIndex = 1:methods
            tic
            for fileIndex = 1:valNum
                data = transpose(tagDeltaMatrix(:, 2, fileIndex));
                switch methodIndex
                    case 1  % weighted with 1/ exp( abs(measurement));
                        [Pbest, result, iters(fileIndex, methodIndex)] = psoSimple4(data, learnRate1, ...
                            learnRate2, inertia, maxIter, dimention, population, ...
                            accuracy, upperBound, lowerBound, convergMaxIter, tagLoc, readerLoc);
                    case 2  % x,y locate seperately, weighted with 1/ exp( abs(measurement))
                        [Pbest, result, iters(fileIndex, methodIndex), result_x2(fileIndex, :), result_y2(fileIndex, :)] = psoSimpleTwin(data, learnRate1, ...
                            learnRate2, inertia, maxIter, dimention, population, ...
                            accuracy, upperBound, lowerBound, convergMaxIter, tagLoc, readerLoc);
                    
                    case 3  %% LSQ optim localization, locate x and y at the same time
                        x0 = transpose([centerPosi(fileIndex, 1) centerPosi(fileIndex, 2) 0]);
                        ydata = zeros((tagPopulation - 2), 1);  % c
                        ydata(1:(tagEachArray - 1)) = tagDeltaMatrix(1:(tagEachArray - 1), 2, fileIndex);
                        ydata(tagEachArray:(tagPopulation - 2)) = tagDeltaMatrix(tagEachArray:(tagPopulation - 2), 2, fileIndex);
                        [x3, resnorm, residual, exitflag, output, lambda, jacobian] = lsqcurvefit(fun3, x0, xdata, ydata, lb, ub);
                        result = [x3(1) x3(2)]; % x should be real location

                    case 4  %% 最小二乘法定位1, locate x and y seperately
                        % x1, x2 is real location and phase bias
                        % x0 is where the fit starts
                        % function output is ideal phase vec calculated using x and xdata
                        % xdata is all the known parameters
                        % ydata is c
                        x0 = transpose([centerPosi(fileIndex, 1) centerPosi(fileIndex, 2) 0]);
                        ydata = tagDeltaMatrix(1:(tagEachArray - 1), 2, fileIndex);
                        [x1, resnorm, residual, exitflag, output, lambda, jacobian] = lsqcurvefit(fun1, x0, xdata, ydata, lb, ub);
                        ydata = tagDeltaMatrix(tagEachArray:(tagPopulation - 2), 2, fileIndex);
                        [x2, resnorm, residual, exitflag, output, lambda, jacobian] = lsqcurvefit(fun2, x0, xdata, ydata, lb, ub);
                        result = [x1(1) x2(2)]; % x should be real location
                end
        %         figure(); plot(Pbest);
                rawResult(fileIndex, :, methodIndex) = result;
            end
            methodIndex
                
            timeVec(methodIndex) = toc;
            % calculate and show accuracy
            [distVec, avaDists(methodIndex, test_index), medDists(methodIndex, test_index),...
                maxDist, variance(methodIndex, test_index)] = ...
                calculateAccuracy(rawResult(:,:,methodIndex), realLoc(:,1:2));
            errorMat(:, methodIndex) = distVec;
        end
        
        if (test_index == 2)
            errorMat2 = errorMat;
        end
        
        timeMatrix(:,test_index) = timeVec;
        [x_error1, y_error1] = AnalyzeXyErr(result_x1, result_y1, realLoc(:, 1:2));
        x_error_vec1(test_index) = x_error1;
        y_error_vec1(test_index) = y_error1;
        [x_error2, y_error2] = AnalyzeXyErr(result_x2, result_y2, realLoc(:, 1:2));
        x_error_vec2(test_index) = x_error2;
        y_error_vec2(test_index) = y_error2;
        % show cummu prab comparison and average localization error of each
        % method
        DrawBarGraph(errorMat);
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
end
%% 
DrawHeightComparison(avaDists, variance);
DrawMethodAccuracyComparison(avaDists, variance, medDists);
DrawFrequencyComparison();
