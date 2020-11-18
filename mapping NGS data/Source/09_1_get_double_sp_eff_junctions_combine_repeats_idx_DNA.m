function get_double_sp_eff_junctions_combine_repeats_idx_DNA_log(idx)

if isstr(idx)
	idx = str2num(idx);
end

out_folder='../Data/Out_double_sp_eff_junctions_idx/';
in_sample_folder='../Data/Out_double_sp_eff_junctions/';

if ~exist(out_folder)
	mkdir(out_folder);
end

files = dir([in_sample_folder '*_DNA_*']);
files={files.name};

files=files(inds);
disp(files)

for i=1:length(files)
    db(i)=load(sprintf('%s%s/Sp_DNA_%d.idx%d.mat',in_sample_folder,files{i},inds(i),idx));
end

double_sp_eff_tbl=db(1).double_tbl_out(:,1);

isoform_1_abundance_mat=zeros(size(double_sp_eff_tbl,1),length(files));
isoform_2_abundance_mat=zeros(size(double_sp_eff_tbl,1),length(files));
isoform_J_abundance_mat=zeros(size(double_sp_eff_tbl,1),length(files));
isoform_B_abundance_mat=zeros(size(double_sp_eff_tbl,1),length(files));

unspliced_abundance_mat=zeros(size(double_sp_eff_tbl,1),length(files));
total_abundance_mat=zeros(size(double_sp_eff_tbl,1),length(files));

for i=1:length(files)
	isoform_1_abundance_mat(:,i)=db(i).double_tbl_out.isoform_1_abundance;
	isoform_2_abundance_mat(:,i)=db(i).double_tbl_out.isoform_2_abundance;
	isoform_J_abundance_mat(:,i)=db(i).double_tbl_out.isoform_J_abundance;
	isoform_B_abundance_mat(:,i)=db(i).double_tbl_out.isoform_B_abundance;
	
	unspliced_abundance_mat(:,i)=db(i).double_tbl_out.unspliced_abundance;
	total_abundance_mat(:,i)=db(i).double_tbl_out.total_abundance;
end

double_sp_eff_tbl.isoform_1_abundance=sum(isoform_1_abundance_mat,2);
double_sp_eff_tbl.isoform_2_abundance=sum(isoform_2_abundance_mat,2);
double_sp_eff_tbl.isoform_J_abundance=sum(isoform_J_abundance_mat,2);
double_sp_eff_tbl.isoform_B_abundance=sum(isoform_B_abundance_mat,2);
double_sp_eff_tbl.unspliced_abundance=sum(unspliced_abundance_mat,2);
double_sp_eff_tbl.total_abundance=sum(total_abundance_mat,2);

double_sp_eff_tbl.isoform_1_sp_eff=(double_sp_eff_tbl.isoform_1_abundance./double_sp_eff_tbl.total_abundance);
double_sp_eff_tbl.isoform_2_sp_eff=(double_sp_eff_tbl.isoform_2_abundance./double_sp_eff_tbl.total_abundance);
double_sp_eff_tbl.isoform_J_sp_eff=(double_sp_eff_tbl.isoform_J_abundance./double_sp_eff_tbl.total_abundance);
double_sp_eff_tbl.isoform_B_sp_eff=(double_sp_eff_tbl.isoform_B_abundance./double_sp_eff_tbl.total_abundance);
double_sp_eff_tbl.unspliced_ratio=(double_sp_eff_tbl.unspliced_abundance./double_sp_eff_tbl.total_abundance);

save(sprintf('%sDNA_junctions_idx%d.mat',out_folder,idx),'double_sp_eff_tbl');



