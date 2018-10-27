function [x_error, y_error] = AnalyzeXyErr(result_x, result_y, realLoc)
x_error1 = median(abs(realLoc(:, 1) - result_x(:,1)));
x_error2 = median(abs(realLoc(:, 2) - result_y(:,2)));
x_error = (x_error1 + x_error2) / 2;
y_error1 = median(abs(realLoc(:, 2) - result_x(:,2)));
y_error2 = median(abs(realLoc(:, 1) - result_y(:,1)));
y_error = (y_error1 + y_error2) / 2;

end