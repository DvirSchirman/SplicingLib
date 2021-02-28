function get_sp_eff_junctions_repeats(currSamplePath, sample_name, currOutPath)


warning('off')

load('double_intron_db.mat')

adaptor_5='TACGTTAAACAATACGAGGCACTTACTCCG';
adaptor_3='CCTGGAGTTCGCTATTCCTATTGTAGTTTT';


filename=sample_name;
load([currSamplePath filename]);

%%
joined_inds=find(cell2mat(cellfun(@(y) ~isempty(strfind(y,'joined')),double_tbl.type,'un',0)));
twointrons_inds=setdiff(1:size(double_tbl,1),joined_inds);
double_tbl=double_tbl(twointrons_inds,:);

[~,inds1,inds2]=intersect(double_tbl.BC,BC_abundance.BC);
double_tbl.abundance(inds1)=BC_abundance.abundance_vec(inds2);
double_tbl.seqs(inds1)=BC_abundance.seq_vec(inds2);

double_tbl.isoform_full=zeros(size(double_tbl,1),1);
double_tbl.isoform_spliced_1=zeros(size(double_tbl,1),1);
double_tbl.isoform_spliced_2=zeros(size(double_tbl,1),1);
double_tbl.isoform_spliced_joined=zeros(size(double_tbl,1),1);
double_tbl.isoform_spliced_both=zeros(size(double_tbl,1),1);

joined_inds=find(cell2mat(cellfun(@(y) ~isempty(strfind(y,'joined')),double_tbl.type,'un',0)));
twointrons_inds=setdiff(1:size(double_tbl,1),joined_inds);
double_tbl=double_tbl(twointrons_inds,:);

double_tbl.seq = cellfun(@(x) [adaptor_5, x, adaptor_3],double_tbl.seq,'un',0);
%%
double_tbl_out=table();
double_tbl_out.BC=double_tbl.BC;


for i=1:size(double_tbl,1)
	if isempty(double_tbl.abundance{i})
			continue
    end
	
	EI5_1_junction=double_tbl.seq{i}(22+double_tbl.SS5_1_inds{i}(1):61+double_tbl.SS5_1_inds{i}(1));
    EI3_1_junction=double_tbl.seq{i}(22+double_tbl.SS3_1_inds{i}(1):61+double_tbl.SS3_1_inds{i}(1));
    EE_1_junction=[double_tbl.seq{i}(22+double_tbl.SS5_1_inds{i}(1):41+double_tbl.SS5_1_inds{i}(1)),...
                    double_tbl.seq{i}(43+double_tbl.SS3_1_inds{i}(end):64+double_tbl.SS3_1_inds{i}(1))];
					
	EI5_2_junction=double_tbl.seq{i}(22+double_tbl.SS5_2_inds{i}(1):61+double_tbl.SS5_2_inds{i}(1));
    EI3_2_junction=double_tbl.seq{i}(22+double_tbl.SS3_2_inds{i}(1):61+double_tbl.SS3_2_inds{i}(1));
    EE_2_junction=[double_tbl.seq{i}(22+double_tbl.SS5_2_inds{i}(1):41+double_tbl.SS5_2_inds{i}(1)),...
                    double_tbl.seq{i}(43+double_tbl.SS3_2_inds{i}(end):64+double_tbl.SS3_2_inds{i}(1))];
					
	EE_J_junction=[double_tbl.seq{i}(22+double_tbl.SS5_1_inds{i}(1):41+double_tbl.SS5_1_inds{i}(1)),...
                    double_tbl.seq{i}(43+double_tbl.SS3_2_inds{i}(end):64+double_tbl.SS3_2_inds{i}(1))];
                
    EE_B_junction=[double_tbl.seq{i}(22+double_tbl.SS5_1_inds{i}(1):41+double_tbl.SS5_1_inds{i}(1)),...
                    double_tbl.seq{i}(43+double_tbl.SS3_1_inds{i}(end):41+double_tbl.SS5_2_inds{i}(1)),...
                    double_tbl.seq{i}(43+double_tbl.SS3_2_inds{i}(end):64+double_tbl.SS3_2_inds{i}(1))];
					
                    
	EI5_1_align_max=swalign(EI5_1_junction,EI5_1_junction);
    EI3_1_align_max=swalign(EI3_1_junction,EI3_1_junction);
    EE_1_align_max=swalign(EE_1_junction,EE_1_junction);
	EI5_2_align_max=swalign(EI5_2_junction,EI5_2_junction);
    EI3_2_align_max=swalign(EI3_2_junction,EI3_2_junction);
    EE_2_align_max=swalign(EE_2_junction,EE_2_junction);
	EE_J_align_max=swalign(EE_J_junction,EE_J_junction);
    EE_B_align_max=swalign(EE_B_junction,EE_B_junction);
	
	EI5_1_alignment = cellfun(@(x) swalign(x,EI5_1_junction,'GapOpen',10,'ExtendGap',1),double_tbl.seqs{i})/EI5_1_align_max;
    EI3_1_alignment = cellfun(@(x) swalign(x,EI3_1_junction,'GapOpen',10,'ExtendGap',1),double_tbl.seqs{i})/EI3_1_align_max;
    EE_1_alignment = cellfun(@(x) swalign(x,EE_1_junction,'GapOpen',10,'ExtendGap',1),double_tbl.seqs{i})/EE_1_align_max;
	EI5_2_alignment = cellfun(@(x) swalign(x,EI5_2_junction,'GapOpen',10,'ExtendGap',1),double_tbl.seqs{i})/EI5_2_align_max;
    EI3_2_alignment = cellfun(@(x) swalign(x,EI3_2_junction,'GapOpen',10,'ExtendGap',1),double_tbl.seqs{i})/EI3_2_align_max;
    EE_2_alignment = cellfun(@(x) swalign(x,EE_2_junction,'GapOpen',10,'ExtendGap',1),double_tbl.seqs{i})/EE_2_align_max;
	EE_J_alignment = cellfun(@(x) swalign(x,EE_J_junction,'GapOpen',10,'ExtendGap',1),double_tbl.seqs{i})/EE_J_align_max;
	EE_B_alignment = cellfun(@(x) swalign(x,EE_B_junction,'GapOpen',10,'ExtendGap',1),double_tbl.seqs{i})/EE_B_align_max;
    
	THR=0.8;
    
    isoform_1_inds = find(EE_1_alignment>THR & EI5_1_alignment<THR & EI3_1_alignment<THR);
    isoform_2_inds = find(EE_2_alignment>THR & EI5_2_alignment<THR & EI3_2_alignment<THR);
    isoform_J_inds = find(EE_J_alignment>THR & EI5_1_alignment<THR & EI3_2_alignment<THR);
    isoform_B_inds = find(EE_B_alignment>THR & EI5_1_alignment<THR & EI3_1_alignment<THR & ...
                        EI5_2_alignment<THR & EI3_2_alignment<THR);
						
	unspliced_inds = find(EE_1_alignment<THR & EI5_1_alignment>THR & EI3_1_alignment>THR & ...
        EE_2_alignment<THR & EI5_2_alignment>THR & EI3_2_alignment>THR );
		
    spliced_inds = [isoform_1_inds; isoform_2_inds; isoform_J_inds; isoform_B_inds];
	other_inds = setdiff(1:length(double_tbl.abundance{i}),[unspliced_inds; spliced_inds]);
	other_inds(double_tbl.abundance{i}(other_inds)<5)=[];
	
	double_tbl_out.isoform_1_abundance(i)=sum(double_tbl.abundance{i}(isoform_1_inds));
	double_tbl_out.isoform_2_abundance(i)=sum(double_tbl.abundance{i}(isoform_2_inds));
	double_tbl_out.isoform_J_abundance(i)=sum(double_tbl.abundance{i}(isoform_J_inds));
	double_tbl_out.isoform_B_abundance(i)=sum(double_tbl.abundance{i}(isoform_B_inds));
	double_tbl_out.unspliced_abundance(i)=sum(double_tbl.abundance{i}(unspliced_inds))+sum(double_tbl.abundance{i}(other_inds));
	
	double_tbl_out.total_abundance(i) = sum([double_tbl_out.isoform_1_abundance(i); double_tbl_out.isoform_2_abundance(i); ...
											double_tbl_out.isoform_J_abundance(i); double_tbl_out.isoform_B_abundance(i); ...
											double_tbl_out.unspliced_abundance(i)]);
    
    a=1;
    
end

sample_str=sample_name(1:end-6);
save([currOutPath sample_str '.mat'],'double_tbl_out');

