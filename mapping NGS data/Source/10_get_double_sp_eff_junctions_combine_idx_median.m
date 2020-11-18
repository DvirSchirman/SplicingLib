close all
clear
clc

in_folder='../Data/Out_double_sp_eff_junctions_idx/';
out_folder='../Data/Splicing_lib/';

prefix = {'DNA';'RNA_log';'RNA_st'};
for k=1:length(prefix)
	files = dir(sprintf('%s/%s*.mat',in_folder,prefix{k}));
	files = {files.name};
	
	for i=1:length(files)
		db(i)=load([in_folder files{i}]);
	end
	double_sp_eff_tbl=db(1).double_sp_eff_tbl(:,1);
	
	isoform_1_mat=zeros(size(double_sp_eff_tbl,1),length(files));
	isoform_2_mat=zeros(size(double_sp_eff_tbl,1),length(files));
	isoform_J_mat=zeros(size(double_sp_eff_tbl,1),length(files));
	isoform_B_mat=zeros(size(double_sp_eff_tbl,1),length(files));
	unspliced_mat=zeros(size(double_sp_eff_tbl,1),length(files));
	total_abundance_mat=zeros(size(double_sp_eff_tbl,1),length(files));
	
	for i=1:length(files)
		isoform_1_mat(:,i)=db(i).double_sp_eff_tbl.isoform_1_sp_eff;
		isoform_2_mat(:,i)=db(i).double_sp_eff_tbl.isoform_2_sp_eff;
		isoform_J_mat(:,i)=db(i).double_sp_eff_tbl.isoform_J_sp_eff;
		isoform_B_mat(:,i)=db(i).double_sp_eff_tbl.isoform_B_sp_eff;
		unspliced_mat(:,i)=db(i).double_sp_eff_tbl.unspliced_ratio;
		total_abundance_mat(:,i)=db(i).double_sp_eff_tbl.total_abundance;
	end
	
	isoform_1(total_abundance_mat<10)=nan;
	isoform_2(total_abundance_mat<10)=nan;
	isoform_J(total_abundance_mat<10)=nan;
	isoform_B(total_abundance_mat<10)=nan;
	unspliced(total_abundance_mat<10)=nan;
	
	isoform_1=nanmedian(isoform_1_mat,2);
	isoform_2=nanmedian(isoform_2_mat,2);
	isoform_J=nanmedian(isoform_J_mat,2);
	isoform_B=nanmedian(isoform_B_mat,2);
	unspliced=nanmedian(unspliced_mat,2);
	
	sum_vec= nansum([isoform_1, isoform_2, isoform_J, isoform_B, unspliced],2);
	
	double_sp_eff_tbl.isoform_1_sp_eff=isoform_1./sum_vec;
	double_sp_eff_tbl.isoform_2_sp_eff=isoform_2./sum_vec;
	double_sp_eff_tbl.isoform_J_sp_eff=isoform_J./sum_vec;
	double_sp_eff_tbl.isoform_B_sp_eff=isoform_B./sum_vec;
	double_sp_eff_tbl.unspliced_ratio=unspliced./sum_vec;
	
	save(sprintf('%s/double_%s_junctions_median.mat',out_folder,prefix{k}),'double_sp_eff_tbl');
end
	