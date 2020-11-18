function get_sp_eff_junctions_combine_repeats_idx_RNA_log(idx)

if isstr(idx)
	idx = str2num(idx);
end

out_folder='../Data/Out_sp_eff_junctions_idx/';
in_sample_folder='../Data/Out_sp_eff_junctions/';

files = dir([in_sample_folder '*_RNA_st*']);
files={files.name};

files=files(inds);
disp(files)

for i=1:length(inds)
    db(i)=load(sprintf('%s%s/Sp_RNA_st_%d.idx%d.mat',in_sample_folder,files{i},inds(i),idx));
end

sp_eff_tbl=db(1).sp_eff_tbl(:,1);
spliced_abundance_mat=zeros(size(sp_eff_tbl,1),length(files));
unspliced_abundance_mat=zeros(size(sp_eff_tbl,1),length(files));
total_abundance_mat=zeros(size(sp_eff_tbl,1),length(files));
for i=1:length(files)
	spliced_abundance_mat(:,i)=db(i).sp_eff_tbl.spliced_abundance;
	unspliced_abundance_mat(:,i)=db(i).sp_eff_tbl.unspliced_abundance;
	total_abundance_mat(:,i)=db(i).sp_eff_tbl.total_abundance;
	cryptic_spliced_abundance_mat(:,i)=db(i).sp_eff_tbl.cryptic_spliced_abundance;
end

sp_eff_tbl.unspliced_abundance=sum(unspliced_abundance_mat,2);
sp_eff_tbl.spliced_abundance=sum(spliced_abundance_mat,2);
sp_eff_tbl.total_abundance=sum(total_abundance_mat,2);
sp_eff_tbl.cryptic_spliced_abundance=sum(cryptic_spliced_abundance_mat,2);

sp_eff_tbl.sp_eff=(sp_eff_tbl.spliced_abundance./sp_eff_tbl.total_abundance);
sp_eff_tbl.cryptic_sp_eff=(sp_eff_tbl.cryptic_spliced_abundance./sp_eff_tbl.total_abundance);

for i=1:size(sp_eff_tbl,1)
	start_inds=[];
	end_inds=[];
	cryptic_abundance_vec=[];
	cryptic_start_inds=[];
	cryptic_end_inds=[];
	cryptic_abundance=[];
	for j=1:length(files)
		start_inds=[start_inds; db(j).sp_eff_tbl.cryptic_start_inds{i}'];
		end_inds=[end_inds; db(j).sp_eff_tbl.cryptic_end_inds{i}'];
		cryptic_abundance_vec=[cryptic_abundance_vec; db(j).sp_eff_tbl.cryptic_inds_abundance{i}'];
	end
	
	start_end_tuple=arrayfun(@(x,y) [num2str(x) ',' num2str(y)],start_inds,end_inds,'un',0);
	tmp = unique(start_end_tuple);
    for j = 1:size(tmp,1)
        start_tmp=str2num(tmp{j}(1:(find(tmp{j}==',',1)-1)));
        end_tmp=str2num(tmp{j}((find(tmp{j}==',',1)+1):end));
        inds=find(start_inds==start_tmp & end_inds==end_tmp);
        cryptic_start_inds(j,1)=start_tmp;
        cryptic_end_inds(j,1)=end_tmp;
        cryptic_abundance(j,1)=sum(cryptic_abundance_vec(inds));
    end

    [cryptic_abundance,inds]=sort(cryptic_abundance,'descend');
    cryptic_start_inds=cryptic_start_inds(inds);
    cryptic_end_inds=cryptic_end_inds(inds);
    
	sp_eff_tbl.cryptic_start_inds{i}=cryptic_start_inds;
	sp_eff_tbl.cryptic_end_inds{i}=cryptic_end_inds;
	sp_eff_tbl.cryptic_inds_abundance{i}=cryptic_abundance;
end

save(sprintf('%sRNA_st_junctions_idx%d.mat',out_folder,idx),'sp_eff_tbl');

	

	