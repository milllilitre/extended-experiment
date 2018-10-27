function DrawHeightComparison(avaDists, variance)
% num_of_lines = size(distMat, 1);
figure(); hold on;
    bar_obj = bar(avaDists(2, :));
    bar_obj.FaceColor = 'flat';
    bar_obj.CData(:,:) = repmat([.5 0 .5], 3, 1);
    errorbar(avaDists(2, :), variance(2, :), 'k', 'lineStyle', 'none');
    
ylabel('Error(m)');
xlabel('Height(m)');
xticks([1 2 3 ]);
xticklabels({'1.3', '1.6', '1.9'});