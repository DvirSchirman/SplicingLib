clear
close all
clc

addpath('../../Data/');
addpath(genpath('../Functions/'));

load('SplicingLib_db.mat')
load('double_site_structure_pairing.mat')
load('splice_site_ids.mat')
load('dG_sliding.mat')
load('GC_sliding.mat')

warning('off')

get_GC = @(x) (sum(x=='G')+sum(x=='C'))/length(x);
get_T = @(x) (sum(x=='T'))/length(x);

splicing_lib_tbl(cellfind(splicing_lib_tbl.type,'double endogenous ss control'),:)=[];
[~,inds1,inds2]=intersect(double_sp_eff_tbl.BC,splicing_lib_tbl.BC);
double_sp_eff_tbl=double_sp_eff_tbl(inds1,:);
splicing_lib_tbl=splicing_lib_tbl(inds2,:);
dG_sliding=dG_sliding(inds2,:);
GC_sliding=GC_sliding(inds2,:);



%%

alternative_splicing_features_tbl=readtable('model_features_all.csv');
alternative_splicing_features_tbl(:,:)=[];
alternative_splicing_features_tbl.Properties.VariableNames{1}='id';

S=42-15;
n=1;
for i=1:size(double_sp_eff_tbl,1)
    alternative_splicing_features_tbl.id{n}=[double_sp_eff_tbl.BC{i} '_1'];
    
    ss5_seq=double_sp_eff_tbl.seq{i}(double_sp_eff_tbl.SS5_1_inds{i}+42);
    ss5_idx=find(strcmp(ss5_map(:,1),ss5_seq));
    if isempty(ss5_idx)
        ss5_idx=size(ss5_map,1);
    end 
    alternative_splicing_features_tbl.ss5_seq(n)=ss5_map{ss5_idx,2};
    
    branch_seq=double_sp_eff_tbl.seq{i}(double_sp_eff_tbl.branch_1_inds{i}+42);
    branch_idx=find(strcmp(branch_map(:,1),branch_seq));
    if isempty(branch_idx)
        branch_idx=size(branch_map,1);
    end
    alternative_splicing_features_tbl.branch_seq(n)=branch_map{branch_idx,2};
    
    ss3_seq=double_sp_eff_tbl.seq{i}(double_sp_eff_tbl.SS3_1_inds{i}+42);
    ss3_idx=find(strcmp(ss3_map(:,1),ss3_seq));
    if isempty(ss3_idx)
        ss3_idx=size(ss3_map,1);
    end 
    alternative_splicing_features_tbl.ss3_seq(n)=ss3_map{ss3_idx,2};
    
    alternative_splicing_features_tbl.splicing_eff_median(n)=double_sp_eff_tbl.isoform_1_sp_eff(i)+double_sp_eff_tbl.isoform_B_sp_eff(i);
    alternative_splicing_features_tbl.intron_len(n)=length(double_sp_eff_tbl.intron_1_seq{i});
    alternative_splicing_features_tbl.intron_GC(n)=get_GC(double_sp_eff_tbl.intron_1_seq{i});
    alternative_splicing_features_tbl.polyT(n)=get_T(double_sp_eff_tbl.intron_1_seq{i}(end-22:end-3));
    
    alternative_splicing_features_tbl.ss5_dG(n)=mean(dG_sliding.dG{i}(S+double_sp_eff_tbl.SS5_1_inds{i}(3:4)));
    alternative_splicing_features_tbl.ss5_GC(n)=mean(GC_sliding.GC{i}(S+double_sp_eff_tbl.SS5_1_inds{i}(3:4)));
    alternative_splicing_features_tbl.ss3_dG(n)=dG_sliding.dG{i}(S+double_sp_eff_tbl.SS3_1_inds{i}(2));
    alternative_splicing_features_tbl.ss3_GC(n)=GC_sliding.GC{i}(S+double_sp_eff_tbl.SS3_1_inds{i}(2));
    alternative_splicing_features_tbl.branch_dG(n)=dG_sliding.dG{i}(S+double_sp_eff_tbl.branch_1_inds{i}(3));
    alternative_splicing_features_tbl.branch_GC(n)=GC_sliding.GC{i}(S+double_sp_eff_tbl.branch_1_inds{i}(3));
    
    alternative_splicing_features_tbl.branch_to_3(n)=double_sp_eff_tbl.SS3_1_inds{i}(end)-double_sp_eff_tbl.branch_1_inds{i}(1)+5;
    
    for j=1:6
        tmp=double_sp_eff_tbl.ss5_1_pairing{i}(j);
        eval(sprintf('alternative_splicing_features_tbl.ss5_pairing_%d(i)=tmp;',j));
    end

    for j=1:3
        tmp=double_sp_eff_tbl.ss3_1_pairing{i}(j);
        eval(sprintf('alternative_splicing_features_tbl.ss3_pairing_%d(i)=tmp;',j))
    end

    for j=1:7
        tmp=double_sp_eff_tbl.branch_1_pairing{i}(j);
        eval(sprintf('alternative_splicing_features_tbl.branch_pairing_%d(i)=tmp;',j))
    end
    
    alternative_splicing_features_tbl.ss5_paired_frac(i)=sum(double_sp_eff_tbl.ss5_1_pairing{i})/6;
    alternative_splicing_features_tbl.branch_paired_frac(i)=sum(double_sp_eff_tbl.branch_1_pairing{i})/7;
    alternative_splicing_features_tbl.ss3_paired_frac(i)=sum(double_sp_eff_tbl.ss3_1_pairing{i})/3;

    alternative_splicing_features_tbl.ss5_stem_len(i)=double_sp_eff_tbl.ss5_1_stem_len(i);
    alternative_splicing_features_tbl.branch_stem_len(i)=double_sp_eff_tbl.branch_1_stem_len(i);
    alternative_splicing_features_tbl.ss5_stem_len(i)=double_sp_eff_tbl.ss3_1_stem_len(i);

    alternative_splicing_features_tbl.ss5_stem_dir(i)=double_sp_eff_tbl.ss5_1_pairing_dir(i);
    alternative_splicing_features_tbl.branch_stem_dir(i)=double_sp_eff_tbl.branch_1_pairing_dir(i);
    alternative_splicing_features_tbl.ss5_stem_dir(i)=double_sp_eff_tbl.ss3_1_pairing_dir(i);
    
    n=n+1;
    
    alternative_splicing_features_tbl.id{n}=[double_sp_eff_tbl.BC{i} '_2'];
    ss5_seq=double_sp_eff_tbl.seq{i}(double_sp_eff_tbl.SS5_2_inds{i}+42);
    ss5_idx=find(strcmp(ss5_map(:,1),ss5_seq));
    if isempty(ss5_idx)
        ss5_idx=size(ss5_map,1);
    end 
    alternative_splicing_features_tbl.ss5_seq(n)=ss5_map{ss5_idx,2};
    
    branch_seq=double_sp_eff_tbl.seq{i}(double_sp_eff_tbl.branch_2_inds{i}+42);
    branch_idx=find(strcmp(branch_map(:,1),branch_seq));
    if isempty(branch_idx)
        branch_idx=size(branch_map,1);
    end
    alternative_splicing_features_tbl.branch_seq(n)=branch_map{branch_idx,2};
    
    ss3_seq=double_sp_eff_tbl.seq{i}(double_sp_eff_tbl.SS3_2_inds{i}+42);
    ss3_idx=find(strcmp(ss3_map(:,1),ss3_seq));
    if isempty(ss3_idx)
        ss3_idx=size(ss3_map,1);
    end 
    alternative_splicing_features_tbl.ss3_seq(n)=ss3_map{ss3_idx,2};
    
    alternative_splicing_features_tbl.splicing_eff_median(n)=double_sp_eff_tbl.isoform_2_sp_eff(i)+double_sp_eff_tbl.isoform_B_sp_eff(i);
    alternative_splicing_features_tbl.intron_len(n)=length(double_sp_eff_tbl.intron_2_seq{i});
    alternative_splicing_features_tbl.intron_GC(n)=get_GC(double_sp_eff_tbl.intron_2_seq{i});
    alternative_splicing_features_tbl.polyT(n)=get_T(double_sp_eff_tbl.intron_2_seq{i}(end-22:end-3));
    
    alternative_splicing_features_tbl.ss5_dG(n)=mean(dG_sliding.dG{i}(S+double_sp_eff_tbl.SS5_2_inds{i}(3:4)));
    alternative_splicing_features_tbl.ss5_GC(n)=mean(GC_sliding.GC{i}(S+double_sp_eff_tbl.SS5_2_inds{i}(3:4)));
    alternative_splicing_features_tbl.ss3_dG(n)=dG_sliding.dG{i}(S+double_sp_eff_tbl.SS3_2_inds{i}(2));
    alternative_splicing_features_tbl.ss3_GC(n)=GC_sliding.GC{i}(S+double_sp_eff_tbl.SS3_2_inds{i}(2));
    alternative_splicing_features_tbl.branch_dG(n)=dG_sliding.dG{i}(S+double_sp_eff_tbl.branch_2_inds{i}(3));
    alternative_splicing_features_tbl.branch_GC(n)=GC_sliding.GC{i}(S+double_sp_eff_tbl.branch_2_inds{i}(3));
    
    alternative_splicing_features_tbl.branch_to_3(n)=double_sp_eff_tbl.SS3_2_inds{i}(end)-double_sp_eff_tbl.branch_2_inds{i}(1)+5;
    
    for j=1:6
        tmp=double_sp_eff_tbl.ss5_2_pairing{i}(j);
        eval(sprintf('alternative_splicing_features_tbl.ss5_pairing_%d(i)=tmp;',j));
    end

    for j=1:3
        tmp=double_sp_eff_tbl.ss3_2_pairing{i}(j);
        eval(sprintf('alternative_splicing_features_tbl.ss3_pairing_%d(i)=tmp;',j))
    end

    for j=1:7
        tmp=double_sp_eff_tbl.branch_2_pairing{i}(j);
        eval(sprintf('alternative_splicing_features_tbl.branch_pairing_%d(i)=tmp;',j))
    end
    
    alternative_splicing_features_tbl.ss5_paired_frac(i)=sum(double_sp_eff_tbl.ss5_2_pairing{i})/6;
    alternative_splicing_features_tbl.branch_paired_frac(i)=sum(double_sp_eff_tbl.branch_2_pairing{i})/7;
    alternative_splicing_features_tbl.ss3_paired_frac(i)=sum(double_sp_eff_tbl.ss3_2_pairing{i})/3;

    alternative_splicing_features_tbl.ss5_stem_len(i)=double_sp_eff_tbl.ss5_2_stem_len(i);
    alternative_splicing_features_tbl.branch_stem_len(i)=double_sp_eff_tbl.branch_2_stem_len(i);
    alternative_splicing_features_tbl.ss5_stem_len(i)=double_sp_eff_tbl.ss3_2_stem_len(i);

    alternative_splicing_features_tbl.ss5_stem_dir(i)=double_sp_eff_tbl.ss5_2_pairing_dir(i);
    alternative_splicing_features_tbl.branch_stem_dir(i)=double_sp_eff_tbl.branch_2_pairing_dir(i);
    alternative_splicing_features_tbl.ss5_stem_dir(i)=double_sp_eff_tbl.ss3_2_pairing_dir(i);
    
    n=n+1;
    
    alternative_splicing_features_tbl.id{n}=[double_sp_eff_tbl.BC{i} '_J'];
    
    ss5_seq=double_sp_eff_tbl.seq{i}(double_sp_eff_tbl.SS5_1_inds{i}+42);
    ss5_idx=find(strcmp(ss5_map(:,1),ss5_seq));
    if isempty(ss5_idx)
        ss5_idx=size(ss5_map,1);
    end 
    alternative_splicing_features_tbl.ss5_seq(n)=ss5_map{ss5_idx,2};
    
    branch_seq=double_sp_eff_tbl.seq{i}(double_sp_eff_tbl.branch_2_inds{i}+42);
    branch_idx=find(strcmp(branch_map(:,1),branch_seq));
    if isempty(branch_idx)
        branch_idx=size(branch_map,1);
    end
    alternative_splicing_features_tbl.branch_seq(n)=branch_map{branch_idx,2};
    
    ss3_seq=double_sp_eff_tbl.seq{i}(double_sp_eff_tbl.SS3_2_inds{i}+42);
    ss3_idx=find(strcmp(ss3_map(:,1),ss3_seq));
    if isempty(ss3_idx)
        ss3_idx=size(ss3_map,1);
    end 
    alternative_splicing_features_tbl.ss3_seq(n)=ss3_map{ss3_idx,2};
    
    alternative_splicing_features_tbl.splicing_eff_median(n)=double_sp_eff_tbl.isoform_J_sp_eff(i);
    alternative_splicing_features_tbl.intron_len(n)=double_sp_eff_tbl.SS3_2_inds{i}(end)-double_sp_eff_tbl.SS5_1_inds{i}(1);
    seq=double_sp_eff_tbl.seq{i}(42+double_sp_eff_tbl.SS5_1_inds{i}(1):42+double_sp_eff_tbl.SS3_2_inds{i}(end));
    alternative_splicing_features_tbl.intron_GC(n)=get_GC(seq);
    alternative_splicing_features_tbl.polyT(n)=get_T(seq(end-22:end-3));
    
    alternative_splicing_features_tbl.ss5_dG(n)=mean(dG_sliding.dG{i}(S+double_sp_eff_tbl.SS5_1_inds{i}(3:4)));
    alternative_splicing_features_tbl.ss5_GC(n)=mean(GC_sliding.GC{i}(S+double_sp_eff_tbl.SS5_1_inds{i}(3:4)));
    alternative_splicing_features_tbl.ss3_dG(n)=dG_sliding.dG{i}(S+double_sp_eff_tbl.SS3_2_inds{i}(2));
    alternative_splicing_features_tbl.ss3_GC(n)=GC_sliding.GC{i}(S+double_sp_eff_tbl.SS3_2_inds{i}(2));
    alternative_splicing_features_tbl.branch_dG(n)=dG_sliding.dG{i}(S+double_sp_eff_tbl.branch_2_inds{i}(3));
    alternative_splicing_features_tbl.branch_GC(n)=GC_sliding.GC{i}(S+double_sp_eff_tbl.branch_2_inds{i}(3));
    
    alternative_splicing_features_tbl.branch_to_3(n)=double_sp_eff_tbl.SS3_2_inds{i}(end)-double_sp_eff_tbl.branch_2_inds{i}(1)+5;
    
    for j=1:6
        tmp=double_sp_eff_tbl.ss5_1_pairing{i}(j);
        eval(sprintf('alternative_splicing_features_tbl.ss5_pairing_%d(i)=tmp;',j));
    end

    for j=1:3
        tmp=double_sp_eff_tbl.ss3_2_pairing{i}(j);
        eval(sprintf('alternative_splicing_features_tbl.ss3_pairing_%d(i)=tmp;',j))
    end

    for j=1:7
        tmp=double_sp_eff_tbl.branch_2_pairing{i}(j);
        eval(sprintf('alternative_splicing_features_tbl.branch_pairing_%d(i)=tmp;',j))
    end
    
    alternative_splicing_features_tbl.ss5_paired_frac(i)=sum(double_sp_eff_tbl.ss5_1_pairing{i})/6;
    alternative_splicing_features_tbl.branch_paired_frac(i)=sum(double_sp_eff_tbl.branch_2_pairing{i})/7;
    alternative_splicing_features_tbl.ss3_paired_frac(i)=sum(double_sp_eff_tbl.ss3_2_pairing{i})/3;

    alternative_splicing_features_tbl.ss5_stem_len(i)=double_sp_eff_tbl.ss5_1_stem_len(i);
    alternative_splicing_features_tbl.branch_stem_len(i)=double_sp_eff_tbl.branch_2_stem_len(i);
    alternative_splicing_features_tbl.ss5_stem_len(i)=double_sp_eff_tbl.ss3_2_stem_len(i);

    alternative_splicing_features_tbl.ss5_stem_dir(i)=double_sp_eff_tbl.ss5_1_pairing_dir(i);
    alternative_splicing_features_tbl.branch_stem_dir(i)=double_sp_eff_tbl.branch_2_pairing_dir(i);
    alternative_splicing_features_tbl.ss5_stem_dir(i)=double_sp_eff_tbl.ss3_2_pairing_dir(i);
     
    n=n+1;
end

alternative_splicing_features_tbl(isnan(alternative_splicing_features_tbl.splicing_eff_median),:)=[];
    
writetable(alternative_splicing_features_tbl,'alternative_sp_model_features_all.csv')
alternative_splicing_features_tbl_binary=alternative_splicing_features_tbl;
alternative_splicing_features_tbl_binary.splicing_eff_median(alternative_splicing_features_tbl_binary.splicing_eff_median>0)=1;
writetable(alternative_splicing_features_tbl_binary,'alternative_sp_model_features_binary_all.csv')

save('alternative_splicing_features.mat','alternative_splicing_features_tbl','alternative_splicing_features_tbl_binary')
%%

all_isoforms_tbl=alternative_splicing_features_tbl(:,1);
all_isoforms_tbl.splicing_eff=alternative_splicing_features_tbl.splicing_eff_median;

alternative_isoforms_tbl=double_sp_eff_tbl(:,cellfind(double_sp_eff_tbl.Properties.VariableNames,'BC'));
alternative_isoforms_tbl.measured_isoform_1=nan*ones(size(alternative_isoforms_tbl,1),1);
alternative_isoforms_tbl.measured_isoform_2=nan*ones(size(alternative_isoforms_tbl,1),1);
alternative_isoforms_tbl.measured_isoform_J=nan*ones(size(alternative_isoforms_tbl,1),1);
alternative_isoforms_tbl.measured_isoform_B=nan*ones(size(alternative_isoforms_tbl,1),1);

for i=1:size(all_isoforms_tbl,1)
    BC=all_isoforms_tbl.id{i}(1:12);
    type=all_isoforms_tbl.id{i}(end);
    ind=cellfind(alternative_isoforms_tbl.BC,BC);
    aa(i)=ind;
    switch type
        case '1'
            alternative_isoforms_tbl.measured_isoform_1(ind)=all_isoforms_tbl.splicing_eff(i);
        case '2'
            alternative_isoforms_tbl.measured_isoform_2(ind)=all_isoforms_tbl.splicing_eff(i);
        case 'J'
            alternative_isoforms_tbl.measured_isoform_J(ind)=all_isoforms_tbl.splicing_eff(i);
    end
end
alternative_isoforms_tbl(isnan(alternative_isoforms_tbl.measured_isoform_J),:)=[];

[~,inds1,inds2]=intersect(alternative_isoforms_tbl.BC,double_sp_eff_tbl.BC);
alternative_isoforms_tbl.measured_isoform_B=double_sp_eff_tbl.isoform_B_sp_eff(inds2);

alternative_isoforms_tbl.state=zeros(size(alternative_isoforms_tbl,1),1);
alternative_isoforms_tbl.state(alternative_isoforms_tbl.measured_isoform_1>0 & ...
    alternative_isoforms_tbl.measured_isoform_2==0 & ...
    alternative_isoforms_tbl.measured_isoform_J==0)=1;
alternative_isoforms_tbl.state(alternative_isoforms_tbl.measured_isoform_1==0 & ...
    alternative_isoforms_tbl.measured_isoform_2>0 & ...
    alternative_isoforms_tbl.measured_isoform_J==0)=2;
alternative_isoforms_tbl.state(alternative_isoforms_tbl.measured_isoform_1>0 & ...
    alternative_isoforms_tbl.measured_isoform_2>0 & ...
    alternative_isoforms_tbl.measured_isoform_J==0)=3;
alternative_isoforms_tbl.state(alternative_isoforms_tbl.measured_isoform_1==0 & ...
    alternative_isoforms_tbl.measured_isoform_2==0 & ...
    alternative_isoforms_tbl.measured_isoform_J>0)=4;
alternative_isoforms_tbl.state(alternative_isoforms_tbl.measured_isoform_1>0 & ...
    alternative_isoforms_tbl.measured_isoform_2==0 & ...
    alternative_isoforms_tbl.measured_isoform_J>0)=5;
alternative_isoforms_tbl.state(alternative_isoforms_tbl.measured_isoform_1==0 & ...
    alternative_isoforms_tbl.measured_isoform_2>0 & ...
    alternative_isoforms_tbl.measured_isoform_J>0)=6;
alternative_isoforms_tbl.state(alternative_isoforms_tbl.measured_isoform_1>0 & ...
    alternative_isoforms_tbl.measured_isoform_2>0 & ...
    alternative_isoforms_tbl.measured_isoform_J>0)=7;

alternative_isoforms_multiclass_tbl=alternative_isoforms_tbl(:,1);
alternative_isoforms_multiclass_tbl.class=alternative_isoforms_tbl.state;

writetable(alternative_isoforms_multiclass_tbl,'alternative_isoforms_multiclass.csv')
