% phase based localization function library
% input and output are row vectors
function output = phaseResolveVec(input)
    % resolve phase, input is a row vector
    length = size(input,2);
    if(length > 0)
        tmpDelta = input(2:length) - input(1:(length - 1));
        tmpT = abs(tmpDelta) > (pi / 2 * 3.25);   %ǰ��ֵ�������Ϊ1������Ϊ0
        tmpT2 = input(1:(length - 1)) > pi;  %ǰֵ����pi��Ϊ1������Ϊ0
        tmpT2 = (tmpT2 - 0.5) * 2; %ǰֵ����pi��Ϊ1������Ϊ-1
        tmpT = tmpT .* tmpT2;   %ǰ�ô�Լpi��������Ϊ1��ǰ��С��pi��������Ϊ-1��ʣ��Ϊ0
        for i = 2:size(tmpT,2)  %�����
            tmpT(i) = tmpT(i - 1) + tmpT(i);
            if(abs(tmpT(i)) == 2)
                tmpT(i) = 0;
            end
        end
        output(1) = input(1);
        output(2:length) = input(2:length) + 2 * pi * tmpT;
    else
        output = 0;
    end
end
