function get_sp_eff_alignment(currSamplePath, sample_name, currOutPath)

load('../Data/double_intron_db.mat')
load('../Data/final_lib.mat')
adaptor_3=seqrcomplement('TAGGAATAGCGAACTCCAGG');

filename=sample_name;
load([currSamplePath filename]);

not_lib_inds=find(cell2mat(cellfun(@(y) ~isempty(strfind(y,'nucleotide composition')),final_lib_tbl.type,'un',0)));
not_lib_inds=[not_lib_inds ;find(cell2mat(cellfun(@(y) ~isempty(strfind(y,'R loops')),final_lib_tbl.type,'un',0)))];
lib_inds = setdiff(1:size(final_lib_tbl,1),not_lib_inds);
lib_tbl=final_lib_tbl(lib_inds,:);

joined_inds=find(cell2mat(cellfun(@(y) ~isempty(strfind(y,'joined')),double_tbl.type,'un',0)));
tmp_tbl=double_tbl(joined_inds,:);
[~,inds1,inds2]=intersect(lib_tbl.BC,tmp_tbl.BC);
lib_tbl.SS5_inds(inds1)=tmp_tbl.SS5_1_inds(inds2);
lib_tbl.SS3_inds(inds1)=tmp_tbl.SS3_1_inds(inds2);
lib_tbl.SS_5(inds1)=cellfun(@(x,y) x(42+y),lib_tbl.seq(inds1),lib_tbl.SS5_inds(inds1),'un',0);
lib_tbl.SS_3(inds1)=cellfun(@(x,y) x(42+y),lib_tbl.seq(inds1),lib_tbl.SS3_inds(inds1),'un',0);

[~,inds1,inds2]=intersect(lib_tbl.BC,BC_abundance.BC);
lib_tbl.abundance(inds1)=BC_abundance.abundance_vec(inds2);
lib_tbl.seqs(inds1)=BC_abundance.seq_vec(inds2);

lib_tbl(cellfind(lib_tbl.type,'Nucleosome occupancy profile'),:)=[];
lib_tbl(find(cellfun(@isempty,lib_tbl.SS5_inds)),:)=[];

%%
lib_tbl.unspliced_abundance=zeros(size(lib_tbl,1),1);
lib_tbl.spliced_abundance=zeros(size(lib_tbl,1),1);
lib_tbl.total_abundance=zeros(size(lib_tbl,1),1);
lib_tbl.cryptic_spliced_abundance=zeros(size(lib_tbl,1),1);
%%

for i=1:size(lib_tbl,1)
	start_inds=[];
    end_inds=[];
    cryptic_abundance=[];
	if isempty(lib_tbl.abundance{i})
			continue
    end
    EI5_junction=lib_tbl.seq{i}(22+lib_tbl.SS5_inds{i}(1):61+lib_tbl.SS5_inds{i}(1));
    EI3_junction=lib_tbl.seq{i}(22+lib_tbl.SS3_inds{i}(1):61+lib_tbl.SS3_inds{i}(1));
    EE_junction=[lib_tbl.seq{i}(22+lib_tbl.SS5_inds{i}(1):41+lib_tbl.SS5_inds{i}(1)),...
                    lib_tbl.seq{i}(43+lib_tbl.SS3_inds{i}(end):64+lib_tbl.SS3_inds{i}(1))];
    
    EI5_align_max=swalign(EI5_junction,EI5_junction);
    EI3_align_max=swalign(EI3_junction,EI3_junction);
    EE_align_max=swalign(EE_junction,EE_junction);
    
    EI5_alignment = cellfun(@(x) swalign(x,EI5_junction,'GapOpen',10,'ExtendGap',1),lib_tbl.seqs{i})/EI5_align_max;
    EI3_alignment = cellfun(@(x) swalign(x,EI3_junction,'GapOpen',10,'ExtendGap',1),lib_tbl.seqs{i})/EI3_align_max;
    EE_alignment = cellfun(@(x) swalign(x,EE_junction,'GapOpen',10,'ExtendGap',1),lib_tbl.seqs{i})/EE_align_max;
    
    THR=0.8;
	spliced_inds = find(EE_alignment>THR & EI5_alignment<THR & EI3_alignment<THR);
    unspliced_inds = find(EE_alignment<THR & EI5_alignment>THR & EI3_alignment>THR);
	other_inds = find(EE_alignment<THR & EI5_alignment<THR & EI3_alignment<THR);
	
	cryptic_inds=[];
    if ~isempty(other_inds)
        ref_seq=lib_tbl.seq{i}(31:end-10);
        ref_alignment=nwalign(ref_seq,ref_seq,'GapOpen',100,'ExtendGap',1);
        [a,b]=cellfun(@(x) nwalign(x,ref_seq,'GapOpen',100,'ExtendGap',1),lib_tbl.seqs{i}(other_inds),'un',0);
        a=cell2mat(a)/ref_alignment;
        intron_inds=cellfun(@(x) find(x(1,:)=='-'),b,'un',0);
        % intron_inds(a>0.7)={[]};
        valid_inds=a<0.7;
		n=1;
        cryptic_start_inds=[];
        cryptic_end_inds=[];
        cryptic_abundance_vec=[];
        for j=1:length(other_inds)
            if ~valid_inds(j)
                continue;
            end
            if sum(diff(intron_inds{j})==1)/length(intron_inds{j})<0.9
                %%% mismatch inds are too sparse
                continue;
            end
            mismatch_inds1=find(b{j}(2,:)==':');
            mismatch_inds2=find(b{j}(2,:)==' ');
            mismatch_inds2=setdiff(mismatch_inds2,intron_inds{j});
            mismatch_inds=[mismatch_inds1, mismatch_inds2];
            if length(mismatch_inds)>6
                %%%too many mismatches
                continue;
            end
			cryptic_start_inds(n,1)=intron_inds{j}(1)-(lib_tbl.SS5_inds{i}(1)+12);
			cryptic_end_inds(n,1)=intron_inds{j}(end)-(lib_tbl.SS3_inds{i}(1)+12);
			cryptic_abundance_vec(n,1)=lib_tbl.abundance{i}(other_inds(j));
            n=n+1;
            cryptic_inds = [cryptic_inds; other_inds(j)];
        end
        start_end_tuple=arrayfun(@(x,y) [num2str(x) ',' num2str(y)],cryptic_start_inds,cryptic_end_inds,'un',0);
        tmp = unique(start_end_tuple);
        for j = 1:size(tmp,1)
            start_tmp=str2num(tmp{j}(1:(find(tmp{j}==',',1)-1)));
            end_tmp=str2num(tmp{j}((find(tmp{j}==',',1)+1):end));
            inds=find(cryptic_start_inds==start_tmp & cryptic_end_inds==end_tmp);
            start_inds(j)=start_tmp;
            end_inds(j)=end_tmp;
            cryptic_abundance(j)=sum(cryptic_abundance_vec(inds));
        end
    end
    
    lib_tbl.unspliced_abundance(i)=sum(lib_tbl.abundance{i}(unspliced_inds));
    lib_tbl.spliced_abundance(i)=sum(lib_tbl.abundance{i}(spliced_inds));
    lib_tbl.cryptic_spliced_abundance(i)=sum(lib_tbl.abundance{i}(cryptic_inds));
    lib_tbl.total_abundance(i)=max([sum(lib_tbl.abundance{i}(lib_tbl.abundance{i}>=5)),...
        lib_tbl.unspliced_abundance(i)+lib_tbl.spliced_abundance(i)+lib_tbl.cryptic_spliced_abundance(i)]);
		
	lib_tbl.cryptic_start_inds{i}=start_inds;
    lib_tbl.cryptic_end_inds{i}=end_inds;
    lib_tbl.cryptic_inds_abundance{i}=cryptic_abundance;    
		
end

% sp_eff_tbl=[lib_tbl.BC, lib_tbl.unspliced_abundance, lib_tbl.predicted_spliced_abundance, lib_tbl.unpredicted_spliced_abundance, lib_tbl.predicted_spliced_seqs,lib_tbl.unpredicted_spliced_seqs];
sp_eff_tbl=[lib_tbl(:,find(strcmp(lib_tbl.Properties.VariableNames,'BC'))),...
			lib_tbl(:,find(strcmp(lib_tbl.Properties.VariableNames,'unspliced_abundance'))),...
			lib_tbl(:,find(strcmp(lib_tbl.Properties.VariableNames,'spliced_abundance'))),...
			lib_tbl(:,find(strcmp(lib_tbl.Properties.VariableNames,'cryptic_spliced_abundance'))),...
			lib_tbl(:,find(strcmp(lib_tbl.Properties.VariableNames,'total_abundance'))),...
            lib_tbl(:,find(strcmp(lib_tbl.Properties.VariableNames,'cryptic_start_inds'))),...
            lib_tbl(:,find(strcmp(lib_tbl.Properties.VariableNames,'cryptic_end_inds'))),...
            lib_tbl(:,find(strcmp(lib_tbl.Properties.VariableNames,'cryptic_inds_abundance'))),...
            
			];
			
			
			
sample_str=sample_name(1:end-6);
save([currOutPath sample_str '.mat'],'sp_eff_tbl');