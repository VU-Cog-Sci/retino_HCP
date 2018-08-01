% load('/Users/martin/Downloads/atlas.mat')
% load('/Users/martin/Downloads/prfresults.mat')

num_v1 = find(startsWith(glasser2016labels,'V1'));
idx_v1 = find(glasser2016==2);

model_fit = 1; % 1 = on all data; 2 = on first half; 3 = on second half
r2_col = find(contains(quants, 'R2'));
sub_avg_r2 = [];
for sub_num = 1:size(subjectids,1)-3 % 3 last are mean
    whole_brain_fit_r2 = allresults(:,r2_col,sub_num,model_fit);
    v1_fit_r2 = whole_brain_fit_r2(idx_v1);
    sub_avg_r2(sub_num,:) = [sub_num,nanmedian(v1_fit_r2)];
end
sort_by_median = sortrows(sub_avg_r2,2,'descend');

num_best = 30;
best_sub_num = sort_by_median(1:num_best,1);
best_sub_label = subjectids(best_sub_num);
