close all
clear
clc

addpath('../Data/');
addpath(genpath('Functions/'));

load('SplicingLib_db.mat')
load('double_intron_db.mat')

tmp=load('double_RNA_log_junctions_median.mat');

[~,inds1,inds2]=intersect(splicing_lib_tbl.BC,tmp.double_sp_eff_tbl.BC);

double_sp_eff_tbl = splicing_lib_tbl(inds1,:);
double_sp_eff_tbl = [double_sp_eff_tbl, tmp.double_sp_eff_tbl(inds2,2:end)];

tmp=load('double_DNA_junctions_median.mat');
DNA_tbl = tmp.double_sp_eff_tbl(inds2,:);

double_sp_eff_tbl.isoform_1_sp_eff(DNA_tbl.isoform_1_sp_eff>0.05)=0;
double_sp_eff_tbl.isoform_2_sp_eff(DNA_tbl.isoform_2_sp_eff>0.05)=0;
double_sp_eff_tbl.isoform_J_sp_eff(DNA_tbl.isoform_J_sp_eff>0.05)=0;
double_sp_eff_tbl.isoform_B_sp_eff(DNA_tbl.isoform_B_sp_eff>0.05)=0;

[~,inds1,inds2] = intersect(double_sp_eff_tbl.BC,double_tbl.BC);

double_sp_eff_tbl.SS5_1_inds(inds1)=double_tbl.SS5_1_inds(inds2);
double_sp_eff_tbl.SS3_1_inds(inds1)=double_tbl.SS3_1_inds(inds2);
double_sp_eff_tbl.SS5_2_inds(inds1)=double_tbl.SS5_2_inds(inds2);
double_sp_eff_tbl.SS3_2_inds(inds1)=double_tbl.SS3_2_inds(inds2);

clear tmp DNA_tbl
save('../Data/double_sp_eff_tbl.mat','double_sp_eff_tbl');

inds = find(cellfun(@isempty, splicing_lib_tbl.SS5_inds));    
splicing_lib_tbl(inds,:)=[];
inds = cellfind(splicing_lib_tbl.type,'R loops - intron');
inds =[inds; cellfind(splicing_lib_tbl.type,'double endogenous - joined ss control')];
inds =[inds; cellfind(splicing_lib_tbl.type,'endogenous - mutated sites')];
inds =[inds; cellfind(splicing_lib_tbl.type,'endogenous pombe - mutated sites')];
inds =[inds; cellfind(splicing_lib_tbl.type,'endogenous ss control')];
inds =[inds; find(cell2mat(cellfun(@(x) ~isempty(strfind(x,'synthetic')),splicing_lib_tbl.type,'un',0)))];
inds =[inds; find(cell2mat(cellfun(@(x) ~isempty(strfind(x,'orthologous')),splicing_lib_tbl.type,'un',0)))];
splicing_lib_tbl(inds,:)=[];

splicing_lib_tbl.intron_seq=cellfun(@(x,y,z) x(y(1)+42:z(end)+42),splicing_lib_tbl.seq,splicing_lib_tbl.SS5_inds,splicing_lib_tbl.SS3_inds,'un',0);

joined_inds = find(cell2mat(cellfun(@(x) ~isempty(strfind(x,'joined')),splicing_lib_tbl.type,'un',0)));
% joined_tbl = splicing_lib_tbl(joined_inds,:);
splicing_lib_tbl(joined_inds,:)=[];

double_sp_eff_tbl.intron_1_seq = cellfun(@(x,y,z) x(42+y(1):42+z(end)),...
    double_sp_eff_tbl.seq,double_sp_eff_tbl.SS5_1_inds,double_sp_eff_tbl.SS3_1_inds,'un',0);
double_sp_eff_tbl.intron_2_seq = cellfun(@(x,y,z) x(42+y(1):42+z(end)),...
    double_sp_eff_tbl.seq,double_sp_eff_tbl.SS5_2_inds,double_sp_eff_tbl.SS3_2_inds,'un',0);
double_sp_eff_tbl.intron_J_seq = cellfun(@(x,y,z) x(42+y(1):42+z(end)),...
    double_sp_eff_tbl.seq,double_sp_eff_tbl.SS5_1_inds,double_sp_eff_tbl.SS3_2_inds,'un',0);

[~,inds1,inds2]=intersect(double_sp_eff_tbl.BC,double_tbl.BC);
double_sp_eff_tbl.branch_1_inds(inds1)=double_tbl.branch_1_inds(inds2);
double_sp_eff_tbl.branch_2_inds(inds1)=double_tbl.branch_2_inds(inds2);
double_sp_eff_tbl.branch_1_seq(inds1)=double_tbl.branch_1_seq(inds2);
double_sp_eff_tbl.branch_2_seq(inds1)=double_tbl.branch_2_seq(inds2);
double_sp_eff_tbl.intron_1_id(inds1)=double_tbl.intron_1_id(inds2);
double_sp_eff_tbl.intron_2_id(inds1)=double_tbl.intron_2_id(inds2);

double_sp_eff_tbl.intron_J_single_sp_eff=nan*ones(size(double_sp_eff_tbl,1),1);
double_sp_eff_tbl.intron_1_single_sp_eff=nan*ones(size(double_sp_eff_tbl,1),1);
double_sp_eff_tbl.intron_2_single_sp_eff=nan*ones(size(double_sp_eff_tbl,1),1);
tic;
for i=1:size(double_sp_eff_tbl,1)
    if ~mod(i,100)
        disp(i)
    end
   
    
    self_align=swalign(double_sp_eff_tbl.intron_1_seq{i},double_sp_eff_tbl.intron_1_seq{i},'GapOpen',100,'ExtendGap',100);
    tmp = cellfun(@(x) swalign(x,double_sp_eff_tbl.intron_1_seq{i},'GapOpen',100,'ExtendGap',100)/self_align,splicing_lib_tbl.intron_seq);
    [a,b]=max(tmp);
    if a>0.8
        double_sp_eff_tbl.intron_1_single_sp_eff(i) = splicing_lib_tbl.splicing_eff_median(b);
    end
   
    self_align=swalign(double_sp_eff_tbl.intron_2_seq{i},double_sp_eff_tbl.intron_2_seq{i},'GapOpen',100,'ExtendGap',100);
    tmp = cellfun(@(x) swalign(x,double_sp_eff_tbl.intron_2_seq{i},'GapOpen',100,'ExtendGap',100)/self_align,splicing_lib_tbl.intron_seq);
    [a,b]=max(tmp);
    if a>0.8
        double_sp_eff_tbl.intron_2_single_sp_eff(i) = splicing_lib_tbl.splicing_eff_median(b);
    end
    
end
toc;
save('../Data/double_sp_eff_tbl.mat','double_sp_eff_tbl')