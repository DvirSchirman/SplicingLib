close all
clear 
clc

addpath('../Data/')

load('FullLib_db.mat')

not_lib_inds=find(cell2mat(cellfun(@(y) ~isempty(strfind(y,'nucleotide composition')),final_lib_tbl.type,'un',0)));
not_lib_inds=[not_lib_inds ;find(cell2mat(cellfun(@(y) ~isempty(strfind(y,'R loops')),final_lib_tbl.type,'un',0)))];
lib_inds = setdiff(1:size(final_lib_tbl,1),not_lib_inds);
lib_tbl=final_lib_tbl(lib_inds,:);

load('RNA_log_junctions_median.mat')
[~,inds1,inds2]=intersect(lib_tbl.BC,sp_eff_tbl.BC);
lib_tbl.splicing_eff_log(inds1)=sp_eff_tbl.sp_eff(inds2);
lib_tbl.cryptic_splicing_eff_log(inds1)=sp_eff_tbl.cryptic_sp_eff(inds2);
lib_tbl.cryptic_start_inds_log(inds1)=sp_eff_tbl.cryptic_start_inds(inds2);
lib_tbl.cryptic_end_inds_log(inds1)=sp_eff_tbl.cryptic_end_inds(inds2);
lib_tbl.cryptic_inds_abundance_log(inds1)=sp_eff_tbl.cryptic_inds_abundance(inds2);
lib_tbl.total_abundance_log(inds1)=sp_eff_tbl.total_abundance(inds2);
lib_tbl.unspliced_ratio_log(inds1)=sp_eff_tbl.unspliced_ratio(inds2);

load('RNA_st_junctions_median.mat')
lib_tbl.splicing_eff_st(inds1)=sp_eff_tbl.sp_eff(inds2);
lib_tbl.cryptic_splicing_eff_st(inds1)=sp_eff_tbl.cryptic_sp_eff(inds2);
lib_tbl.cryptic_start_inds_st(inds1)=sp_eff_tbl.cryptic_start_inds(inds2);
lib_tbl.cryptic_end_inds_st(inds1)=sp_eff_tbl.cryptic_end_inds(inds2);
lib_tbl.cryptic_inds_abundance_st(inds1)=sp_eff_tbl.cryptic_inds_abundance(inds2);
lib_tbl.total_abundance_st(inds1)=sp_eff_tbl.total_abundance(inds2);
lib_tbl.unspliced_ratio_st(inds1)=sp_eff_tbl.unspliced_ratio(inds2);

load('DNA_junctions_median.mat')
[~,inds1,inds2]=intersect(lib_tbl.BC,sp_eff_tbl.BC);
lib_tbl.splicing_eff_DNA(inds1)=sp_eff_tbl.sp_eff(inds2);
lib_tbl.cryptic_splicing_eff_DNA(inds1)=sp_eff_tbl.cryptic_sp_eff(inds2);
lib_tbl.cryptic_start_inds_DNA(inds1)=sp_eff_tbl.cryptic_start_inds(inds2);
lib_tbl.cryptic_end_inds_DNA(inds1)=sp_eff_tbl.cryptic_end_inds(inds2);
lib_tbl.cryptic_inds_abundance_DNA(inds1)=sp_eff_tbl.cryptic_inds_abundance(inds2);
lib_tbl.expression_log(inds1)=log10((lib_tbl.total_abundance_log(inds1)/nansum(lib_tbl.total_abundance_log))./(sp_eff_tbl.total_abundance(inds2)/nansum(sp_eff_tbl.total_abundance)));
lib_tbl.expression_st(inds1)=log10((lib_tbl.total_abundance_st(inds1)/nansum(lib_tbl.total_abundance_st))./(sp_eff_tbl.total_abundance(inds2)/nansum(sp_eff_tbl.total_abundance)));

nucleosome_inds=find(cell2mat(cellfun(@(y) ~isempty(strfind(y,'Nucleosome')),lib_tbl.type,'un',0)));
nucleosome_lib_tbl=lib_tbl(nucleosome_inds,:);
splicing_lib_tbl=lib_tbl(setdiff(1:size(lib_tbl,1),nucleosome_inds),:);

inds = find(splicing_lib_tbl.splicing_eff_DNA>0.05);
splicing_lib_tbl.splicing_eff_st(inds)=0;
splicing_lib_tbl.splicing_eff_log(inds)=0;

inds = find(splicing_lib_tbl.cryptic_splicing_eff_DNA>0.05 & splicing_lib_tbl.cryptic_splicing_eff_log>0);
for i= 1:length(inds)
    ind = inds(i);
    inds_tuple = arrayfun(@(x,y) [num2str(x) ',' num2str(y)],splicing_lib_tbl.cryptic_start_inds_log{ind},splicing_lib_tbl.cryptic_end_inds_log{ind},'un',0);
    inds_tuple_DNA = arrayfun(@(x,y) [num2str(x) ',' num2str(y)],splicing_lib_tbl.cryptic_start_inds_DNA{ind},splicing_lib_tbl.cryptic_end_inds_DNA{ind},'un',0);
    [~,inds_rm,~]=intersect(inds_tuple,inds_tuple_DNA);
    rm_abundance=sum(splicing_lib_tbl.cryptic_inds_abundance_log{ind}(inds_rm));
    splicing_lib_tbl.cryptic_splicing_eff_log(ind)=splicing_lib_tbl.cryptic_splicing_eff_log(ind)*(1-rm_abundance/sum(splicing_lib_tbl.cryptic_inds_abundance_log{ind}));
    if length(inds_rm)==length(splicing_lib_tbl.cryptic_end_inds_log{ind})
        splicing_lib_tbl.cryptic_end_inds_log{ind}=[];
        splicing_lib_tbl.cryptic_start_inds_log{ind}=[];
        splicing_lib_tbl.cryptic_inds_abundance_log{ind}=[];
    else
        splicing_lib_tbl.cryptic_end_inds_log{ind}(inds_rm)=[];
        splicing_lib_tbl.cryptic_start_inds_log{ind}(inds_rm)=[];
        splicing_lib_tbl.cryptic_inds_abundance_log{ind}(inds_rm)=[];
    end
end

inds = find(splicing_lib_tbl.cryptic_splicing_eff_DNA>0.05 & splicing_lib_tbl.cryptic_splicing_eff_st>0);
for i= 1:length(inds)
    ind = inds(i);
    inds_tuple = arrayfun(@(x,y) [num2str(x) ',' num2str(y)],splicing_lib_tbl.cryptic_start_inds_st{ind},splicing_lib_tbl.cryptic_end_inds_st{ind},'un',0);
    inds_tuple_DNA = arrayfun(@(x,y) [num2str(x) ',' num2str(y)],splicing_lib_tbl.cryptic_start_inds_DNA{ind},splicing_lib_tbl.cryptic_end_inds_DNA{ind},'un',0);
    [~,inds_rm,~]=intersect(inds_tuple,inds_tuple_DNA);
    rm_abundance=sum(splicing_lib_tbl.cryptic_inds_abundance_st{ind}(inds_rm));
    splicing_lib_tbl.cryptic_splicing_eff_st(ind)=splicing_lib_tbl.cryptic_splicing_eff_st(ind)*(1-rm_abundance/sum(splicing_lib_tbl.cryptic_inds_abundance_st{ind}));
    if length(inds_rm)==length(splicing_lib_tbl.cryptic_end_inds_st{ind})
        splicing_lib_tbl.cryptic_end_inds_st{ind}=[];
        splicing_lib_tbl.cryptic_start_inds_st{ind}=[];
        splicing_lib_tbl.cryptic_inds_abundance_st{ind}=[];
    else
        splicing_lib_tbl.cryptic_end_inds_st{ind}(inds_rm)=[];
        splicing_lib_tbl.cryptic_start_inds_st{ind}(inds_rm)=[];
        splicing_lib_tbl.cryptic_inds_abundance_st{ind}(inds_rm)=[];
    end
end

splicing_lib_tbl.splicing_eff_median=splicing_lib_tbl.splicing_eff_log;

inds=find(cell2mat(cellfun(@(x) ~isempty(strfind(x,'DNA')),splicing_lib_tbl.Properties.VariableNames,'un',0)));
splicing_lib_tbl(:,inds)=[];

splicing_lib_tbl.splicing_eff_median(splicing_lib_tbl.total_abundance_log<10)=nan;

%%

%%
inds = find(cellfun(@(x) ~isempty(strfind(x,'joined')),splicing_lib_tbl.type));
inds = [inds; find(cellfun(@(x) ~isempty(strfind(x,'pombe')),splicing_lib_tbl.type))];
inds = [inds; find(cellfun(@(x) ~isempty(strfind(x,'intron control')),splicing_lib_tbl.type))];
inds = setdiff(inds,cellfind(splicing_lib_tbl.type,'double endogenous intron control'));
splicing_lib_tbl(inds,:)=[];

save('../Data/SplicingLib_db.mat','splicing_lib_tbl')




