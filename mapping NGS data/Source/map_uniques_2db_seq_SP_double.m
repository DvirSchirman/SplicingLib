function map_uniques_2db_R(currSamplePath, sample_name, currOutPath)
% clear
% currSamplePath='../Sequencing_data/RNA/Splicing_lib/Out_derep_untrimmed/Sp_RNA_st_3/';
% sample_name='Sp_RNA_st_3.idx1.fasta';
% currOutPath='../Sequencing_data/RNA/Splicing_lib/Out_mapped_seq_double/Sp_RNA_st_3/';

db=load('double_intron_db.mat'); 
% db=load('Nucleotide_comp_db.mat'); 
% db=load('R_loops_db.mat'); 

% samplePath='../Sequencing_data/preliminary_sequencing/Splicing_lib/Out_derep/';
% outPath='../Sequencing_data/preliminary_sequencing/Splicing_lib/Out_mapped/';
progPath=[currOutPath 'progress/'];

% if ~isdir(outPath)
	% mkdir(outPath);
% end
if ~isdir(progPath)
	mkdir(progPath);
end


%%%% read BC & abundance from fasta
filename=sample_name;
S=fastaread([currSamplePath filename]);
tmp=struct2cell(S);
seq=tmp(2,:);
BC=cellfun(@(x) x(1:12),seq,'un',0);
tmp=strsplit([S.Header],{';','='});
tmp=tmp(3:3:length(tmp));
abundance=cellfun(@str2num,tmp);

fid=fopen([progPath sample_name '.txt'],'w');

%%%% map direct matches to db
n=1;
BC_abundance.BC{n}='XXXXX';
BC_abundance.abundance_vec{n}=1;
BC_abundance.seq_vec{n}='a';
n=n+1;
for i=1:length(db.double_tbl.BC)
    inds=find(strcmp(BC,db.double_tbl.BC{i}));
	if ~isempty(inds)
		BC_abundance.BC{n}=db.double_tbl.BC{i};
		for j=1:length(inds)
			BC_abundance.abundance_vec{n}(j,1)=abundance(inds(j));
			BC_abundance.seq_vec{n}(j,1)=seq(inds(j));
		end
        n=n+1;
    end
	fprintf(fid,'direct match: %d\n',i);
end

%%%% map minimal edit distance
BC_len=12;
count=0;
for i=1:length(BC)
	if ismember(BC{i},BC_abundance.BC)
        continue;
    end
	% len=length(BC{i});
    
	[m_f,ind_f]=sort(EditDistance2(BC{i},cell2mat(db.double_tbl.BC)));
	% BC_tmp=db.double_tbl.BC(BC_len==len);
	BC_min=db.double_tbl.BC(ind_f(1));
    
	
	% [m_f,ind_f]=sort([m{1}(1),m{2}(1),m{3}(1)]);
	if m_f(1)<3 && m_f(1)<m_f(2)
		count=count+1;
		BC_tmp=BC_min;
		if ismember(BC_tmp,BC_abundance.BC)
			tmp_ind=find(strcmp(BC_tmp,BC_abundance.BC));
			BC_abundance.abundance_vec{tmp_ind}(end+1)=abundance(i);
			BC_abundance.seq_vec{tmp_ind}(end+1)=seq(i);
			%BC_abundance_tab(ind).LIG1=BC_abundance.LIG1(ind);
		else
			BC_abundance.BC(end+1)=BC_tmp;
			BC_abundance.abundance_vec{end+1}=abundance(i);
			BC_abundance.seq_vec{end+1}=seq(i);
		end
	end
	fprintf(fid,'distance match: %d\n',i);
end
fclose(fid);
BC_abundance.BC(1)=[];
BC_abundance.abundance_vec(1)=[];
BC_abundance.seq_vec(1)=[];
sample_str=sample_name(1:end-4);
save([currOutPath sample_str '.mat'],'BC_abundance');




