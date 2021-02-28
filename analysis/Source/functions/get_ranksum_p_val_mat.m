function p_mat = get_ranksum_p_val_mat(samples_cell)

p_mat = nan*ones(length(samples_cell));
for i = 1:length(samples_cell)
    for j= i+1:length(samples_cell)
%         [~, p_mat(i,j)]=ttest2(samples_cell{i},samples_cell{j});
        if sum(~isnan(samples_cell{i}))<2 | sum(~isnan(samples_cell{j}))<2
            p_mat(i,j)=nan;
        else
            p_mat(i,j)=ranksum(samples_cell{i},samples_cell{j});
        end
    end
end