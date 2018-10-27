% phase based localization function library
% input and output are row vectors
% only used as a means to regulate phase into a rising slope
% 允许测量值非单调递增
function output = phaseRegulateVec(input)
    % resolve phase, input is a row vector
    length = size(input,2);
    if(length > 0)
        tmpDelta = input(2:length) - input(1:(length - 1));
        tmpT = (tmpDelta < (- pi * 0.5)) - (tmpDelta > (pi * 1.5));
        for i = 2:size(tmpT,2)  %求积分
            tmpT(i) = tmpT(i - 1) + tmpT(i);
        end
        output(1) = input(1);
        output(2:length) = input(2:length) + 2 * pi * tmpT;
    else
        output = 0;
    end
end
