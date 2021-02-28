clear
close all
clc

addpath('../Data/');
addpath(genpath('Functions/'));

Figures_str = '../Figures/Figure5/';
if ~exist(Figures_str)
    mkdir(Figures_str)
end

load('SplicingLib_db.mat')
load('unintended_splicing_tbl.mat')

get_GC = @(x) (sum(x=='G')+sum(x=='C'))/length(x);

ortholog_inds = find(cell2mat(cellfun(@(y) (strcmp(y,'orthologous')),splicing_lib_tbl.type,'un',0)));
ortholog_inds = [ortholog_inds ; find(cell2mat(cellfun(@(y) (strcmp(y,'endogenous')),splicing_lib_tbl.type,'un',0)))];
ortholog_tbl = splicing_lib_tbl(ortholog_inds,:);
ortholog_tbl(cellfun(@length, ortholog_tbl.branch_pos_from_3)~=1,:)=[];
ortholog_tbl.branch_pos_from_3=cell2mat(ortholog_tbl.branch_pos_from_3);
ortholog_tbl(cellfind(ortholog_tbl.orth_organism,'agos'),:)=[];


ortholog_tbl.GC=cellfun(@(x) get_GC(x), ortholog_tbl.seq);
ortholog_tbl.intron_seq=cellfun(@(x,y,z) x(y(1)+42:z(end)+42),ortholog_tbl.seq,ortholog_tbl.SS5_inds,ortholog_tbl.SS3_inds,'un',0);
ortholog_tbl.intron_GC=cellfun(@(x) get_GC(x), ortholog_tbl.intron_seq);

ortholog_SE_tbl=table();

ortholog_SE_tbl.gene=unique(ortholog_tbl.intron_cerevisiae_gene);

species_list={'cerevisiae';'cgla';'ncas';'ndai';'kafr';'tbla';'ecym';'zrou';'tdel';'kthe';'klac'};

for i=1:length(species_list)
    ortholog_SE_tbl.(species_list{i})=-0.5*ones(size(ortholog_SE_tbl,1),1);
end

ortholog_dSE_tbl=ortholog_SE_tbl(:,[1,3:end]);
ortholog_dSE_tbl{:,2:end}=-1.5*ones(size(ortholog_dSE_tbl{:,2:end}));

for i=1:size(ortholog_tbl,1)
    if strcmp(ortholog_tbl.type{i},'orthologous')
        ind=find(strcmp(ortholog_tbl.intron_cerevisiae_gene,ortholog_tbl.intron_cerevisiae_gene{i}) & strcmp(ortholog_tbl.type,'endogenous'));
        ortholog_tbl.diff_sp_eff(i)=ortholog_tbl.splicing_eff_median(i)-nanmean(ortholog_tbl.splicing_eff_median(ind));
    end
end

for i=1:size(ortholog_tbl,1)
    gene = ortholog_tbl.intron_cerevisiae_gene{i};
    species = ortholog_tbl.orth_organism{i};
    if isempty(species)
        species = 'cerevisiae';
    end
    ind = find(cell2mat(cellfun(@(x) strcmp(x,gene),ortholog_SE_tbl.gene,'un',0)));
    eval(sprintf('ortholog_SE_tbl.%s(ind)=ortholog_tbl.splicing_eff_median(i);',species));
    if strcmp(ortholog_tbl.type{i},'orthologous')
        eval(sprintf('ortholog_dSE_tbl.%s(ind)=ortholog_tbl.diff_sp_eff(i);',species));
    end
end

%%
species_list={'S.cerevisiae';'C.glabrata';'N.castellii';'N.dairenensis';'K.africana';'T.blattae';'E.cymbalariae';'Z.rouxii';'T.delbrueckii';'K.thermotolerans';'K.lactis'};
SE_mat=ortholog_SE_tbl{:,2:end};
SE_mat(SE_mat==-0.5)=nan;

%% Figure 5A

figure('units','centimeters','outerposition',[2 2 27 16])
myviolinplot_mean(mat2cell(SE_mat,140,ones(1,size(SE_mat,2))))
ax = gca;
for i=1:length(species_list)
    label_str{i} = sprintf('%s [N=%d]',species_list{i},sum(~isnan(SE_mat(:,i))));
end
ax.XTickLabel = label_str;
ax.XTickLabelRotation=-45;
ax.FontSize = 25;

export_fig(sprintf('%sA - species_SE.png',Figures_str),'-png','-r100','-transparent');

%% Figure 5B

d_SE_mat=ortholog_dSE_tbl{:,2:end};
d_SE_mat(d_SE_mat==-1.5)=nan;
figure('units','centimeters','outerposition',[2 2 27 13.5])
hold on
CT=cbrewer('qual','Set1',size(d_SE_mat,2)+2);
CT(8,:)=[];
for i=1:size(d_SE_mat,2)
    bar([zeros(1,i) sum(d_SE_mat(:,i)>0)./sum(~isnan(d_SE_mat(:,i)))*100 zeros(1,size(d_SE_mat,2)-i)],'FaceColor',CT(i+1,:));
end
ylim([0 100])
ax=gca;
ax.XTick=1:size(d_SE_mat,2)+1;
ax.XTickLabel=[{''},species_list(2:end)'];
ax.XTickLabelRotation=-45;
ax.FontSize = 25;
xlim([0 size(d_SE_mat,2)+2])

export_fig(sprintf('%sB - fraction_of_pos_dSE.png',Figures_str),'-png','-r100','-transparent');

%% Figure 5C

u2af1_species = {'zrou','tdel','ecym','kthe'}; 
[~,u2af1_inds,~]=intersect(ortholog_SE_tbl.Properties.VariableNames,u2af1_species);
no_u2af1_inds=setdiff(2:size(ortholog_SE_tbl,2),u2af1_inds);
no_u2af1_inds(no_u2af1_inds==7)=[]; % remove blattae because it has mutated U2AF1

figure
u2af1_SE=SE_mat(:,u2af1_inds-1);
no_u2af1_SE=SE_mat(:,no_u2af1_inds-1);
myviolinplot_mean([{u2af1_SE(:)},{no_u2af1_SE(:)}])
[~,p]=ttest2(u2af1_SE(:),no_u2af1_SE(:));
ax = gca;
ax.XTickLabel = {'U2AF1(+)','U2AF1(-)'};
ax.FontSize = 22;
ylabel('Splicing efficiency','fontsize',18)

export_fig(sprintf('%sC - u2af1.png',Figures_str),'-png','-r100','-transparent');

%% Figure 5D

load('ss3_avoidance_W30_combined.mat')

p_window_downstream_vec(p_window_downstream_vec==0)=1e-5;
p_window_upstream_vec(p_window_upstream_vec==0)=1e-5;

organisms={'S.cerevisiae';'C.glabrata';'N.castellii';'N.dairenensis';'K.africana';'T.blattae';'E.cymbalariae';'Z.rouxii';'T.delbrueckii';'K.thermotolerans';'K.lactis'};
inds=[1,2,4,6,3,7,10,8,9,11,5];
p_window_downstream_vec=p_window_downstream_vec(inds);
p_window_upstream_vec=p_window_upstream_vec(inds);

organisms_plus_inds=[6,7,8,9,10];
organisms_minus_inds=setdiff(1:length(organisms),organisms_plus_inds);
organisms_plus_inds(1)=[];

figure('units','centimeters','outerposition',[2 2 27 27])
subplot(2,1,1)
hold on
bar(-log10(p_window_upstream_vec))
plot(0:length(organisms)+1,-log10(5e-2/(1*length(p_window_downstream_vec))*ones(1,length(p_window_downstream_vec)+2))','r--','linewidth',2.5)
ax = gca;
ax.XTick=1:length(organisms);
ax.XTickLabel=[];
ax.FontSize = 25;

subplot(2,1,2)
hold on
bar(-log10(p_window_downstream_vec))
plot(0:length(organisms)+1,-log10(5e-2/(1*length(p_window_downstream_vec))*ones(1,length(p_window_downstream_vec)+2))','r--','linewidth',2.5)
ax = gca;
ax.XTick=1:length(organisms);
ax.XTickLabel=organisms;
ax.XTickLabelRotation=360-45;
ax.FontSize = 25;

export_fig(sprintf('%sD - ss3_avoidance_p.png',Figures_str),'-png','-r100','-transparent');