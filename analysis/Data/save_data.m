close all
clear
clc

load('SplicingLib_db.mat')

get_GC =@(x) (sum(x=='C')+sum(x=='G'))./length(x);
get_Y = @(x) (sum(x=='T')+sum(x=='C'))/length(x);
get_T = @(x) (sum(x=='T'))/length(x);

inds=find(cellfun(@isempty,splicing_lib_tbl.intron_len));
splicing_lib_tbl(inds,:)=[];
inds=find(cellfun(@length,splicing_lib_tbl.branch_pos_from_3)~=1);
splicing_lib_tbl(inds,:)=[];
inds=find(cellfun(@length,splicing_lib_tbl.W_tail_len)~=1);
splicing_lib_tbl(inds,:)=[];

splicing_lib_tbl.intron_seq=cellfun(@(x,y,z) x(y(1)+42:z(end)+42),splicing_lib_tbl.seq,splicing_lib_tbl.SS5_inds,splicing_lib_tbl.SS3_inds,'un',0);
splicing_lib_tbl.intron_len = cell2mat(splicing_lib_tbl.intron_len);
splicing_lib_tbl.branch_pos_from_3 = cell2mat(splicing_lib_tbl.branch_pos_from_3);
splicing_lib_tbl.W_tail_len = cell2mat(splicing_lib_tbl.W_tail_len);
splicing_lib_tbl.polyT=cellfun(@(x) get_T(x(end-19:end)), splicing_lib_tbl.intron_seq);

splicing_lib_tbl(:,[6,7,8,10,19,20,21,27,30,33,38,39,40,41,42,44,52:58, 60, 61])=[];
tmp = splicing_lib_tbl.Properties.VariableNames;
for i=1:length(tmp)
    if length(tmp{i})>3
        if strcmp(tmp{i}(end-3:end),'_log')
            tmp{i} = tmp{i}(1:end-4);
        end
    end
end
splicing_lib_tbl.Properties.VariableNames = tmp;

inds = find(cellfun(@(x) ~isempty(x),splicing_lib_tbl.SS5_inds));
splicing_lib_tbl.SS5_inds(inds) = cellfun(@(x) x(1),splicing_lib_tbl.SS5_inds(inds),'un',0);
inds = find(cellfun(@(x) ~isempty(x),splicing_lib_tbl.SS3_inds));
splicing_lib_tbl.SS3_inds(inds) = cellfun(@(x) x(1),splicing_lib_tbl.SS3_inds(inds),'un',0);
inds = find(cellfun(@(x) ~isempty(x),splicing_lib_tbl.branch_inds));
splicing_lib_tbl.branch_inds(inds) = cellfun(@(x) x(1),splicing_lib_tbl.branch_inds(inds),'un',0);
inds = find(cellfun(@(x) ~isempty(x),splicing_lib_tbl.cryptic_start_inds));

splicing_lib_tbl(:,[13,15,31:33])=[];
writetable(splicing_lib_tbl,'SplicingLibData.csv')