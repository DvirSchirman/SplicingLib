close all
clear
clc

load('double_sp_eff_tbl.mat')

get_GC =@(x) (sum(x=='C')+sum(x=='G'))./length(x);
get_Y = @(x) (sum(x=='T')+sum(x=='C'))/length(x);
get_T = @(x) (sum(x=='T'))/length(x);

double_sp_eff_tbl = removevars(double_sp_eff_tbl,{'DRIP','Gamma','Tx','source_gene','SS_5','branch_seq','SS_3','intron_len','branch_pos_from_3',...
    'W_tail_len','site_hairpin','background','len','NL_profile','seqN','prob','mutated_sites','mutations_per_site','SS5_inds','branch_inds',...
    'SS3_inds','ss3_stem_percent','ss3_dG','NL_profile_type','ss3_stem_len','branch_stem_percent','branch_dG','branch_stem_len','ss5_stem_percent',...
    'ss5_dG','ss5_stem_len','intron_cerevisiae_gene','orth_organism','cost_Chen','cost_Lynch','predicted_NL_profile','dG','MTTR',...
    'splicing_eff_log','splicing_eff_median','splicing_eff_st','total_abundance_log','total_abundance_st','unspliced_ratio_st','expression_st',...
    'unspliced_ratio_log','nucleotide_composition','GC','expression_log'});

double_sp_eff_tbl(:,cellfun(@(x) ~isempty(strfind(x,'cryptic')),double_sp_eff_tbl.Properties.VariableNames)) = [];

double_sp_eff_tbl.SS5_1_inds = cellfun(@(x) x(1),double_sp_eff_tbl.SS5_1_inds);
double_sp_eff_tbl.SS5_2_inds = cellfun(@(x) x(1),double_sp_eff_tbl.SS5_2_inds);
double_sp_eff_tbl.SS3_1_inds = cellfun(@(x) x(1),double_sp_eff_tbl.SS3_1_inds);
double_sp_eff_tbl.SS3_2_inds = cellfun(@(x) x(1),double_sp_eff_tbl.SS3_2_inds);
double_sp_eff_tbl.branch_1_inds = cellfun(@(x) x(1),double_sp_eff_tbl.branch_1_inds);
double_sp_eff_tbl.branch_2_inds = cellfun(@(x) x(1),double_sp_eff_tbl.branch_2_inds);

double_sp_eff_tbl.SS5_1_seq = cellfun(@(x,y) x(42+y:42+y+5),double_sp_eff_tbl.seq,mat2cell(double_sp_eff_tbl.SS5_1_inds,ones(size(double_sp_eff_tbl,1),1)),'un',0);
double_sp_eff_tbl.SS3_1_seq = cellfun(@(x,y) x(42+y:42+y+2),double_sp_eff_tbl.seq,mat2cell(double_sp_eff_tbl.SS3_1_inds,ones(size(double_sp_eff_tbl,1),1)),'un',0);
double_sp_eff_tbl.SS5_2_seq = cellfun(@(x,y) x(42+y:42+y+5),double_sp_eff_tbl.seq,mat2cell(double_sp_eff_tbl.SS5_2_inds,ones(size(double_sp_eff_tbl,1),1)),'un',0);
double_sp_eff_tbl.SS3_2_seq = cellfun(@(x,y) x(42+y:42+y+2),double_sp_eff_tbl.seq,mat2cell(double_sp_eff_tbl.SS3_2_inds,ones(size(double_sp_eff_tbl,1),1)),'un',0);
double_sp_eff_tbl = movevars(double_sp_eff_tbl,{'SS5_1_seq'},'Before','branch_1_seq');
double_sp_eff_tbl = movevars(double_sp_eff_tbl,{'SS3_1_seq'},'After','branch_1_seq');
double_sp_eff_tbl = movevars(double_sp_eff_tbl,{'SS5_2_seq'},'Before','branch_2_seq');
double_sp_eff_tbl = movevars(double_sp_eff_tbl,{'SS3_2_seq'},'After','branch_2_seq');


double_sp_eff_tbl = movevars(double_sp_eff_tbl,{'BC','isoform_1_sp_eff','isoform_2_sp_eff','isoform_J_sp_eff','isoform_B_sp_eff'},'After',1);
double_sp_eff_tbl = movevars(double_sp_eff_tbl,{'branch_1_inds'},'After','SS5_1_inds');
double_sp_eff_tbl = movevars(double_sp_eff_tbl,{'branch_2_inds'},'After','SS5_2_inds');

writetable(double_sp_eff_tbl,'SplicingLibData_DoubleIntrons_noIndex.xlsx')