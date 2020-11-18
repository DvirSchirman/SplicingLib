clear
close all
clc

addpath('../../Data/');
addpath(genpath('../Functions/'));

load('SplicingLib_db.mat')

get_GC = @(x) (sum(x=='G')+sum(x=='C'))/length(x);
get_T = @(x) (sum(x=='T'))/length(x);


%%
load('double_intron_db.mat')
inds=find(cell2mat(cellfun(@(x) ~isempty(strfind(x,'joined')),double_tbl.type,'un',0)));
double_tbl=double_tbl(inds,:);
double_tbl(cellfind(double_tbl.type,'double endogenous - joined ss control'),:)=[];
double_tbl(:,cellfind(double_tbl.Properties.VariableNames,'SS5_inds'))=[];
double_tbl(:,cellfind(double_tbl.Properties.VariableNames,'SS3_inds'))=[];

double_tbl.Properties.VariableNames{cellfind(double_tbl.Properties.VariableNames,'SS5_1_inds')}='SS5_inds';
double_tbl.Properties.VariableNames{cellfind(double_tbl.Properties.VariableNames,'SS3_1_inds')}='SS3_inds';
double_tbl.SS_5=cellfun(@(x,y) x(y+12),double_tbl.seq,double_tbl.SS5_inds,'un',0);
double_tbl.SS_3=cellfun(@(x,y) x(y+12),double_tbl.seq,double_tbl.SS3_inds,'un',0);
double_tbl.intron_len=cellfun(@(x,y) y(end)-x(1),double_tbl.SS5_inds,double_tbl.SS3_inds);

branch_inds=cellfun(@(x) strfind(x,'TAAC')-3-12,double_tbl.seq,'un',0);
double_tbl(find(cellfun(@isempty,branch_inds)),:)=[];
branch_inds(find(cellfun(@isempty,branch_inds)))=[];
branch_from_3=cellfun(@(x,y) y(end)-x,branch_inds,double_tbl.SS3_inds,'un',0);
for i=1:length(branch_inds)
    if length(branch_inds{i})>1
        [a,b]=min(abs(branch_from_3{i}-15));
        branch_inds{i}=branch_inds{i}(b);
        branch_from_3{i}=branch_from_3{i}(b);
    end
end
double_tbl(find(cellfun(@(x) x<1, branch_inds)),:)=[];
branch_inds(find(cellfun(@(x) x<1, branch_inds)))=[];

branch_seqs=cellfun(@(x,y) x(y+12:y+12+6),double_tbl.seq,branch_inds,'un',0);

double_tbl.branch_inds=branch_inds;
double_tbl.branch_seq=branch_seqs;
%%

inds = find(cellfun(@isempty, splicing_lib_tbl.SS5_inds)); 
inds = unique([inds; find(cellfun(@isempty, splicing_lib_tbl.SS_5))]);
splicing_lib_tbl(inds,:)=[];

ctrl_inds = find(cell2mat(cellfun(@(x) ~isempty(strfind(x,'control')),splicing_lib_tbl.type,'un',0)));
ctrl_inds = [ctrl_inds; cellfind(splicing_lib_tbl.orth_organism,'agos')];
ctrl_tbl=splicing_lib_tbl(ctrl_inds,:);

splicing_lib_tbl=splicing_lib_tbl(setdiff(1:size(splicing_lib_tbl,1),ctrl_inds),:);

splicing_lib_tbl.intron_seq=cellfun(@(x,y,z) x(y(1)+42:z(end)+42),splicing_lib_tbl.seq,splicing_lib_tbl.SS5_inds,splicing_lib_tbl.SS3_inds,'un',0);
splicing_lib_tbl.intron_GC=cellfun(@(x) get_GC(x), splicing_lib_tbl.intron_seq);
% splicing_lib_tbl.intron_GC(splicing_lib_tbl.intron_GC==0)=nan;
splicing_lib_tbl.polyT=cellfun(@(x) get_T(x(end-22:end-3)), splicing_lib_tbl.intron_seq);

features_tbl=splicing_lib_tbl(:,[cellfind(splicing_lib_tbl.Properties.VariableNames,'BC'),cellfind(splicing_lib_tbl.Properties.VariableNames,'splicing_eff_median')]);

load('site_dG_GC_30.mat')
load('site_structure_pairing.mat')

[~,inds1,inds2]=intersect(features_tbl.BC,sliding_dG_GC_tbl.BC);
features_tbl=features_tbl(inds1,:);
splicing_lib_tbl=splicing_lib_tbl(inds1,:);
sliding_dG_GC_tbl=sliding_dG_GC_tbl(inds2,:);
[~,inds1,inds2]=intersect(features_tbl.BC,site_pairing_tbl.BC);
site_pairing_tbl=site_pairing_tbl(inds2,:);

%%

ss5_seqs=splicing_lib_tbl.SS_5;
ss5_u=tabulate(ss5_seqs);
[~,inds]=sort(cell2mat(ss5_u(:,2)),'descend');
ss5_u=ss5_u(inds,:);
tmp_ind=find(cell2mat(ss5_u(:,2))<10,1);
ss5_map=ss5_u(1:tmp_ind,1);
% ss5_map(:,2)=mat2cell(randperm(size(ss5_map,1))',ones(size(ss5_map,1),1));
ss5_map(:,2)=mat2cell(flipud((1:size(ss5_map,1))'),ones(size(ss5_map,1),1));
ss5_map(tmp_ind+1:size(ss5_u,1),1)=ss5_u(tmp_ind+1:size(ss5_u,1),1);
ss5_map(tmp_ind+1:size(ss5_u,1),2)=ss5_map(tmp_ind,2);
features_tbl.ss5_seq=cellfun(@(x) ss5_map(find(strcmp(ss5_map(:,1),x)),2),ss5_seqs);


tmp_branch_inds=find(cell2mat(cellfun(@(y) (~isempty(strfind(y,'NN'))),splicing_lib_tbl.branch_seq,'un',0)));
tmp_branch_tbl=splicing_lib_tbl(tmp_branch_inds,:);
tmp_branch_seqs=cellfun(@(x,y) x(y(1)+40:y(1)+46),tmp_branch_tbl.seq,tmp_branch_tbl.branch_inds,'un',0);
splicing_lib_tbl.branch_seq(tmp_branch_inds)=tmp_branch_seqs;

branch_seqs=splicing_lib_tbl.branch_seq;
branch_seqs(cellfun(@isempty,branch_seqs))={'ZZZZZZZ'};
branch_u=tabulate(branch_seqs);
[~,inds]=sort(cell2mat(branch_u(:,2)),'descend');
branch_u=branch_u(inds,:);
tmp_ind=find(cell2mat(branch_u(:,2))<10,1);
branch_map=branch_u(1:tmp_ind,1);
% branch_map(:,2)=mat2cell(randperm(size(branch_map,1))',ones(size(branch_map,1),1));
branch_map(:,2)=mat2cell(flipud((1:size(branch_map,1))'),ones(size(branch_map,1),1));
branch_map(tmp_ind+1:size(branch_u,1),1)=branch_u(tmp_ind+1:size(branch_u,1),1);
branch_map(tmp_ind+1:size(branch_u,1),2)=branch_map(tmp_ind,2);
features_tbl.branch_seq=cellfun(@(x) branch_map(find(strcmp(branch_map(:,1),x)),2),branch_seqs);


ss3_seqs=splicing_lib_tbl.SS_3;
ss3_u=tabulate(ss3_seqs);
[~,inds]=sort(cell2mat(ss3_u(:,2)),'descend');
ss3_u=ss3_u(inds,:);
tmp_ind=find(cell2mat(ss3_u(:,2))<10,1);
ss3_map=ss3_u(1:tmp_ind,1);
% ss3_map(:,2)=mat2cell(randperm(size(ss3_map,1))',ones(size(ss3_map,1),1));
ss3_map(:,2)=mat2cell(flipud((1:size(ss3_map,1))'),ones(size(ss3_map,1),1));
ss3_map(tmp_ind+1:size(ss3_u,1),1)=ss3_u(tmp_ind+1:size(ss3_u,1),1);
ss3_map(tmp_ind+1:size(ss3_u,1),2)=ss3_map(tmp_ind,2);
features_tbl.ss3_seq=cellfun(@(x) ss3_map(find(strcmp(ss3_map(:,1),x)),2),ss3_seqs);


features_tbl.intron_GC=splicing_lib_tbl.intron_GC;
features_tbl.polyT=splicing_lib_tbl.polyT;
features_tbl.intron_len=splicing_lib_tbl.intron_len;
features_tbl.branch_to_3=splicing_lib_tbl.branch_pos_from_3;
features_tbl.branch_to_3(cellfun(@length,features_tbl.branch_to_3)>1)={27};
features_tbl.ss5_dG=sliding_dG_GC_tbl.ss5_dG;
features_tbl.ss5_GC=sliding_dG_GC_tbl.ss5_GC;
features_tbl.branch_dG=sliding_dG_GC_tbl.branch_dG;
features_tbl.branch_GC=sliding_dG_GC_tbl.branch_GC;
features_tbl.ss3_dG=sliding_dG_GC_tbl.ss3_dG;
features_tbl.ss3_GC=sliding_dG_GC_tbl.ss3_GC;

for i=1:6
    features_tbl(:,end+1)=cellfun(@(x) x(i),site_pairing_tbl.ss5_pairing,'un',0);
    features_tbl.Properties.VariableNames{end}=sprintf('ss5_pairing_%d',i);
end
for i=1:7
    features_tbl(:,end+1)=cellfun(@(x) x(i),site_pairing_tbl.branch_pairing,'un',0);
    features_tbl.Properties.VariableNames{end}=sprintf('branch_pairing_%d',i);
end
for i=1:3
    features_tbl(:,end+1)=cellfun(@(x) x(i),site_pairing_tbl.ss3_pairing,'un',0);
    features_tbl.Properties.VariableNames{end}=sprintf('ss3_pairing_%d',i);
end
features_tbl.ss5_paired_frac=cellfun(@sum,site_pairing_tbl.ss5_pairing)/6;
features_tbl.branch_paired_frac=cellfun(@sum,site_pairing_tbl.branch_pairing)/7;
features_tbl.ss3_paired_frac=cellfun(@sum,site_pairing_tbl.ss3_pairing)/3;

features_tbl.ss5_stem_len=site_pairing_tbl.ss5_stem_len;
features_tbl.branch_stem_len=site_pairing_tbl.branch_stem_len;
features_tbl.ss3_stem_len=site_pairing_tbl.ss3_stem_len;
features_tbl.ss5_stem_dir=site_pairing_tbl.ss5_pairing_dir;
features_tbl.ss5_stem_dir(features_tbl.ss5_stem_dir<0)=2;
features_tbl.branch_stem_dir=site_pairing_tbl.branch_pairing_dir;
features_tbl.branch_stem_dir(features_tbl.branch_stem_dir<0)=2;
features_tbl.ss3_stem_dir=site_pairing_tbl.ss3_pairing_dir;
features_tbl.ss3_stem_dir(features_tbl.ss3_stem_dir<0)=2;

%%

features_tbl = features_tbl(randperm(size(features_tbl,1)),:);
features_tbl=features_tbl(isfinite(features_tbl.splicing_eff_median),:);


save('features_tbl_all.mat','features_tbl');
writetable(features_tbl,'model_features_all.csv')


%%
sites_freqs=zeros(size(features_tbl,1),3);

SS_5=['GTATGT';'GTACGT';'GTATGA';'GTAAGT';'GTATGC';'GCATGT'];
SS_5_freq=[198; 24; 10; 9; 9 ; 4];
SS_5_freq=SS_5_freq/sum(SS_5_freq);
ss5_map(:,3)={0};
for i=1:size(SS_5,1)
    ind = cellfind(ss5_map(:,1),SS_5(i,:));
    ss5_map{ind,3}=SS_5_freq(i);
    inds = find(cell2mat(features_tbl.ss5_seq)==ss5_map{ind,2});
    sites_freqs(inds,1)=SS_5_freq(i);
end

SS_3=['TAG';'CAG';'AAG'];
SS_3_freq=[137;116;7];
SS_3_freq=SS_3_freq/sum(SS_3_freq);
ss3_map(:,3)={0};
for i=1:size(SS_3,1)
    ind = cellfind(ss3_map(:,1),SS_3(i,:));
    ss3_map{ind,3}=SS_3_freq(i);
    inds = find(cell2mat(features_tbl.ss3_seq)==ss3_map{ind,2});
    sites_freqs(inds,3)=SS_3_freq(i);
end

branch=['NNCTAAC';'NNCTAAT';'NNTTAAC';'TACTAAC'];
branch_freq=[253-226,1,4,226];
branch_freq=branch_freq/sum(branch_freq(1:4));

branch_map_5 = cellfun(@(x) x(3:end),branch_map(:,1),'un',0);
for i=1:size(branch,1)
    if i==1
        ind1 = cellfind(branch_map_5,branch(i,3:end));
        ind2 = cellfind(branch_map(:,1),branch(4,:));
        ind = setdiff(ind1,ind2);
    elseif i==4
        ind = cellfind(branch_map(:,1),branch(4,:));
    else
        ind = cellfind(branch_map_5,branch(i,3:end));
    end
        
    branch_map(ind,3)={branch_freq(i)};
    for j=1:length(ind)
        inds = find(cell2mat(features_tbl.branch_seq)==branch_map{ind(j),2});
        sites_freqs(inds,2)=branch_freq(i);
    end
end
inds = find(cellfun(@isempty,branch_map(:,3)));
branch_map(inds,3)={0};

% writetable(array2table(sites_freqs),'../gradient boosting/sites_frequencies.csv')
tbl=cell2table(ss5_map);
tbl.Properties.VariableNames={'seq','id','freq'};
writetable(tbl,'ss5_map.csv')
tbl=cell2table(ss3_map);
tbl.Properties.VariableNames={'seq','id','freq'};
writetable(tbl,'ss3_map.csv')
tbl=cell2table(branch_map);
tbl.Properties.VariableNames={'seq','id','freq'};
writetable(tbl,'branch_map.csv')
save('splice_site_ids.mat','ss5_map','branch_map','ss3_map')




