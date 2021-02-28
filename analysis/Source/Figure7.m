clear
close all
clc

addpath('../Data/');
addpath(genpath('Functions/'));

Figures_str = '../Figures/Figure7/';
if ~exist(Figures_str)
    mkdir(Figures_str)
end

load('SplicingLib_db.mat')
load('double_sp_eff_tbl.mat')

%% Figure 7B
double_sp_eff_tbl.isoform_1_sp_eff(double_sp_eff_tbl.isoform_1_sp_eff<0.05)=0;
double_sp_eff_tbl.isoform_2_sp_eff(double_sp_eff_tbl.isoform_2_sp_eff<0.05)=0;
double_sp_eff_tbl.isoform_J_sp_eff(double_sp_eff_tbl.isoform_J_sp_eff<0.05)=0;
double_sp_eff_tbl.isoform_B_sp_eff(double_sp_eff_tbl.isoform_B_sp_eff<0.05)=0;

figure('units','centimeters','outerposition',[2 2 32 18])
myviolinplot_mean([{double_sp_eff_tbl.isoform_1_sp_eff(double_sp_eff_tbl.isoform_1_sp_eff>0)},...
                {double_sp_eff_tbl.isoform_2_sp_eff(double_sp_eff_tbl.isoform_2_sp_eff>0)},...
                {double_sp_eff_tbl.isoform_J_sp_eff(double_sp_eff_tbl.isoform_J_sp_eff>0)},...
                {double_sp_eff_tbl.isoform_B_sp_eff(double_sp_eff_tbl.isoform_B_sp_eff>0)}])
            
xlim([0.5 4.5])            
ax = gca;
ax.FontSize=35;

export_fig(sprintf('%sB - isoforms_SE.png',Figures_str),'-png','-r100','-transparent');

%% Figure 7C

introns = {};
for i=1:12
    if i<11
        introns{i}=sprintf('cerevisiae_%d',i);
    else
        introns{i}=sprintf('pombe_%d',i-2);
    end
end
for i=13:24
    if i<23
        introns{i}=sprintf('intron_ctrl_cerevisiae_%d',i-12);
    else
        introns{i}=sprintf('intron_ctrl_pombe_%d',i-2-12);
    end
end


for i=1:length(introns)
    inds1 = cellfind(double_sp_eff_tbl.intron_1_id,introns{i});
    introns_isoform1_spliced_ratio(i) = sum(double_sp_eff_tbl.isoform_1_sp_eff(inds1)>0)/length(inds1);
    introns_isoform1_mean_sp_eff(i) = nanmean(double_sp_eff_tbl.isoform_1_sp_eff(inds1));
    
    inds2 = cellfind(double_sp_eff_tbl.intron_2_id,introns{i});
    introns_isoform2_spliced_ratio(i) = sum(double_sp_eff_tbl.isoform_2_sp_eff(inds2)>0)/length(inds2);
    introns_isoform2_mean_sp_eff(i) = nanmean(double_sp_eff_tbl.isoform_2_sp_eff(inds2));
    
    isoform_mat_1=[double_sp_eff_tbl.isoform_1_sp_eff(inds1),double_sp_eff_tbl.isoform_2_sp_eff(inds1),double_sp_eff_tbl.isoform_J_sp_eff(inds1)];
    isoform_mat_2=[double_sp_eff_tbl.isoform_1_sp_eff(inds2),double_sp_eff_tbl.isoform_2_sp_eff(inds2),double_sp_eff_tbl.isoform_J_sp_eff(inds2)];
    
    M=(length(inds1)+length(inds2));
    unspliced(i) = sum(sum([isoform_mat_1; isoform_mat_2]>0,2)==0)/M;
    self_only(i) = (sum((isoform_mat_1(:,1)>0).*(sum(isoform_mat_1>0,2)==1)) + ...
                sum((isoform_mat_2(:,2)>0).*(sum(isoform_mat_2>0,2)==1)))/M;
    self_plus_others(i) = (sum((isoform_mat_1(:,1)>0).*(sum(isoform_mat_1>0,2)>1)) + ...
                sum((isoform_mat_2(:,2)>0).*(sum(isoform_mat_2>0,2)>1)))/M;
    others_only(i) = (sum((isoform_mat_1(:,1)==0).*(sum(isoform_mat_1>0,2)>0)) + ...
                sum((isoform_mat_2(:,2)==0).*(sum(isoform_mat_2>0,2)>0)))/M;
            
    isoform_B(i) = (sum(double_sp_eff_tbl.isoform_B_sp_eff(inds1)>0)+sum(double_sp_eff_tbl.isoform_B_sp_eff(inds2)>0))/M;

end


figure
scatter(introns_isoform1_spliced_ratio,introns_isoform2_spliced_ratio,100,'filled','MarkerFaceAlpha',0.7)
[r,p]=nancorr(introns_isoform1_spliced_ratio',introns_isoform2_spliced_ratio','Spearman')
xlabel('Intron 1 - proportion spliced','fontsize',20)
ylabel('Intron 2 - proportion spliced','fontsize',20)
ax = gca;
ax.FontSize = 22;

export_fig(sprintf('%sC - double_intron_proportion_scatter.png',Figures_str),'-png','-r100','-transparent');

%% Figure 7D

inds=1:10;

introns_tmp=introns(inds);
% introns_sp_eff_tmp=introns_sp_eff(inds);
introns_sp_eff_tmp=self_only(inds)+self_plus_others(inds);
 
class_mat=nan*ones(length(introns_tmp));
[introns_sp_eff_tmp, inds]=sort(introns_sp_eff_tmp);
introns_tmp=introns_tmp(inds);
isoform_mat=[double_sp_eff_tbl.isoform_1_sp_eff,double_sp_eff_tbl.isoform_2_sp_eff,double_sp_eff_tbl.isoform_J_sp_eff];
isoform_mat=isoform_mat>0;
double_sp_eff_tbl.class = isoform_mat*[1;2;4];
for i=1:length(introns_tmp)
    for j=1:length(introns_tmp)
        inds1 = cellfind(double_sp_eff_tbl.intron_1_id,introns_tmp{i});
        inds2 = cellfind(double_sp_eff_tbl.intron_2_id,introns_tmp{j});
        inds=intersect(inds1,inds2);
        if length(inds)>1
            tmp=sum(isoform_mat(inds,:)>0)>0;
            class_mat(i,j)=tmp(1)+2*tmp(2)+4*tmp(3);
        elseif length(inds)==1
            class_mat(i,j)=double_sp_eff_tbl.class(inds);
        end
    end
end

figure
class_mat(find(class_mat==3 | class_mat>=5))=6; 
class_mat(find(class_mat<5 & class_mat>0))=2; 
imagesc(flipud(class_mat))
CT = cbrewer('qual','Set1',9);
CT = CT([9,2,1],:);
colormap(CT)
ax = gca;
ax.YTick = 1:length(introns_tmp);
ax.XTick = 1:length(introns_tmp);
ax.YTickLabel=length(introns_tmp):-1:1;
ax.FontSize = 25;

ylabel('Intron 1','fontsize',16)  
xlabel('Intron 2','fontsize',16)

export_fig(sprintf('%sD - double_intron_isoform_num_intron_rank.png',Figures_str),'-png','-r100','-transparent');

%% Figure 6E
% load('Gradient boosting (Figures 5 + 6E)/alternative_splicing_features.mat')
% 
% features_mat=alternative_splicing_features_tbl{:,3:end};
% alternative_splicing_features_multiclass=alternative_splicing_features_tbl(1:3:end,1);
% alternative_splicing_features_multiclass.Properties.VariableNames{1}='BC';
% alternative_splicing_features_multiclass.BC = cellfun(@(x) x(1:12),alternative_splicing_features_multiclass.BC,'un',0);
% 
% features_mat = [features_mat(1:3:end,:), features_mat(2:3:end,:) ,features_mat(3:3:end,:)];
% features_mat = features_mat(:,[1:76, 80,81]);
% 
% alternative_splicing_features_multiclass{:,2:size(features_mat,2)+1} = features_mat;
% alternative_splicing_features_multiclass.Properties.VariableNames(2:39)=...
%     cellfun(@(x) [x '_1'],alternative_splicing_features_tbl.Properties.VariableNames(3:end),'un',0);
% alternative_splicing_features_multiclass.Properties.VariableNames(40:77)=...
%     cellfun(@(x) [x '_2'],alternative_splicing_features_tbl.Properties.VariableNames(3:end),'un',0);
% alternative_splicing_features_multiclass.Properties.VariableNames(78:79)={'intron_GC_J','intron_len_J'};
% 
% clear alternative_splicing_features_tbl alternative_splicing_features_tbl_binary
% 
% shap_values=zeros(size(features_mat,1),size(features_mat,2),8);
% for i=1:8
%     tmp=readtable(sprintf('Gradient boosting (Figures 5 + 6E)/alternative_isoform_shap_values_multiclass_%d.csv',i-1));
%     shap_values(:,:,i)=table2array(tmp);
% end
% 
% inds=[0,1,2,4,3,5,6,7]+1;
% shap_values=shap_values(:,:,inds);
% shap_values_mean = squeeze(mean(abs(shap_values),1));
% 
% [a,b]=find(shap_values_mean>3*std(shap_values(:)));
% inds=unique(a);
% 
% ref_shap_mat=table2array(readtable('Gradient boosting (Figures 5 + 6E)/shap_values_isoforms/shap_reference.csv'));
% ref_shap_mat=[ref_shap_mat, ref_shap_mat(:,[2,3])];
% ref_features_mat=table2array(readtable('Gradient boosting (Figures 5 + 6E)/shap_values_isoforms/features_reference.csv'));
% ref_features_mat=[ref_features_mat, ref_features_mat(:,[2,3])];
% ref_corr=diag(nancorr(ref_shap_mat,ref_features_mat));
% ref_corr_sign=ref_corr>0;
% 
% for i=1:8
%     tmp_shap=table2array(readtable(sprintf('Gradient boosting (Figures 5 + 6E)/shap_values_isoforms/isoform_shap_values_multiclass_%d.csv',i-1)));
%     tmp_features=table2array(readtable(sprintf('Gradient boosting (Figures 5 + 6E)/shap_values_isoforms/isoform_feature_values_multiclass_%d.csv',i-1)));
%     tmp_corr=diag(nancorr(tmp_shap,tmp_features));
%     tmp_corr_sign=tmp_corr>0;
%     sign_mat(:,i)=-2*xor(tmp_corr_sign,ref_corr_sign)+1;
% end
% 
% inds1=[0,1,2,4,3,5,6,7]+1;
% sign_mat=sign_mat(:,inds1);
% 
% figure('units','centimeters','outerposition',[2 2 38 18])
% imagesc(shap_values_mean(inds,:).*sign_mat)
% CT=(cbrewer('div','RdBu',256));
% % CT=flipud(CT);
% colormap(CT)
% c=colorbar;
% ax = gca;
% ax.CLim = [-0.4 0.4];
% ax.YTick = (1:length(inds));
% str=alternative_splicing_features_multiclass.Properties.VariableNames(inds+1);
% for i=1:length(str)
%     str{i}(str{i}=='_')=' ';
% end
% ax.YTickLabel = str;
% ax.XTick=1:8;
% ax.XTickLabel={'Unspliced','Isoform 1','Isoform 2','Isoform J','Isoform 1+2','Isoform 1+J','Isoform 2+J','Isoform 1+2+J'};
% ax.FontSize = 16;
% ylabel(c,'Effect-sign * |mean(shap values)|','fontsize',16);    
% 
% export_fig(sprintf('%sE - alternative_isoforms_SHAP.png',Figures_str),'-png','-r100','-transparent');
% 
% 
% 
