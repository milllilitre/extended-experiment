function DrawBarGraph(distMat)
% num_of_lines = size(distMat, 1);
num_of_methods = size(distMat, 2);
ava_vec = zeros(1, num_of_methods);
var_vec = zeros(1, num_of_methods);
figure(); hold on;
for method_index = 1:num_of_methods
    ava_vec(method_index) = mean(distMat(:, method_index));
    var_vec(method_index) = var(distMat(:, method_index));
    bar(ava_vec);
    errorbar(ava_vec, var_vec, 'k', 'LineStyle', 'none');
end
% legend('show');
ylabel('Error(m)');
xticks([1 2.1 3 4]);
xticklabels({'PSO+weight', '2PSO+weight', 'LSQ', '2LSQ'});