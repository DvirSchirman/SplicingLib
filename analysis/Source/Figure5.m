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

get_GC =@(x) (sum(x=='C')+sum(x=='G'))./length(x);

inds = find(cellfun(@isempty, splicing_lib_tbl.SS5_inds));    
splicing_lib_tbl(inds,:)=[];
inds = find(cellfun(@isempty, splicing_lib_tbl.SS3_inds));    
splicing_lib_tbl(inds,:)=[];

ctrl_inds = find(cell2mat(cellfun(@(x) ~isempty(strfind(x,'control')),splicing_lib_tbl.type,'un',0)));
ctrl_inds = [ctrl_inds; cellfind(splicing_lib_tbl.orth_organism,'agos')];
ctrl_tbl=splicing_lib_tbl(ctrl_inds,:);

splicing_lib_tbl=splicing_lib_tbl(setdiff(1:size(splicing_lib_tbl,1),ctrl_inds),:);

splicing_lib_tbl.GC=cellfun(@(x) get_GC(x),splicing_lib_tbl.seq);
splicing_lib_tbl.intron_seq=cellfun(@(x,y,z) x(y(1)+42:z(end)+42),splicing_lib_tbl.seq,splicing_lib_tbl.SS5_inds,splicing_lib_tbl.SS3_inds,'un',0);
splicing_lib_tbl.intron_GC=cellfun(@(x) get_GC(x), splicing_lib_tbl.intron_seq);

load('site_dG_GC_30.mat')

splicing_lib_tbl.ss3_dG=nan*ones(size(splicing_lib_tbl,1),1);
[~,inds1,inds2]=intersect(splicing_lib_tbl.BC,sliding_dG_GC_tbl.BC);
splicing_lib_tbl.ss3_dG(inds1)=sliding_dG_GC_tbl.ss3_dG(inds2);

inds_11=find(splicing_lib_tbl.splicing_eff_median>0 & splicing_lib_tbl.cryptic_splicing_eff_log>0);
inds_01=find(splicing_lib_tbl.splicing_eff_median==0 & splicing_lib_tbl.cryptic_splicing_eff_log>0);
inds_10=find(splicing_lib_tbl.splicing_eff_median>0 & splicing_lib_tbl.cryptic_splicing_eff_log==0);
inds_00=find(splicing_lib_tbl.splicing_eff_median==00 & splicing_lib_tbl.cryptic_splicing_eff_log==0);

%% Figure 5A
start_inds = cell2mat(splicing_lib_tbl.cryptic_start_inds_log);
end_inds = cell2mat(splicing_lib_tbl.cryptic_end_inds_log);
abundances = cell2mat(splicing_lib_tbl.cryptic_inds_abundance_log);

% inds_tuple = arrayfun(@(x,y) [num2str(x) ',' num2str(y)],start_inds,end_inds,'un',0);
start_inds_u = unique(start_inds);
start_inds_abundances=zeros(length(start_inds_u),1);
for i=1:length(start_inds_u)
    inds=find(start_inds==start_inds_u(i));
    start_inds_abundances(i)=sum(abundances(inds));
end

end_inds_u = unique(end_inds);
end_inds_abundances=zeros(length(end_inds_u),1);
for i=1:length(end_inds_u)
    inds=find(end_inds==end_inds_u(i));
    end_inds_abundances(i)=sum(abundances(inds));
end

W=3;
low_bin=round((min([start_inds_u;end_inds_u])+W/2)/W,0)*W;
high_bin=round((max([start_inds_u;end_inds_u])-W/2)/W,0)*W;
bins=low_bin:W:high_bin;
start_bins_abundances=zeros(length(bins),1);
end_bins_abundances=zeros(length(bins),1);
for i=1:length(bins)
    if i==1
        inds=find(start_inds_u<bins(i)+W/2);
    elseif i==length(bins)
        inds=find(start_inds_u>=bins(i)-W/2);
    else
        inds=find(start_inds_u>=bins(i)-W/2 & start_inds_u<bins(i)+W/2);
    end
    if ~isempty(inds)
        start_bins_abundances(i)=sum(start_inds_abundances(inds));
    end
    
    if i==1
        inds=find(end_inds_u<bins(i)+W/2);
    elseif i==length(bins)
        inds=find(end_inds_u>=bins(i)-W/2);
    else
        inds=find(end_inds_u>=bins(i)-W/2 & end_inds_u<bins(i)+W/2);
    end
    if ~isempty(inds)
        end_bins_abundances(i)=sum(end_inds_abundances(inds));
    end
end

figure('units','centimeters','outerposition',[2 2 32 13.5])
subplot(1,2,1)
bar(bins,start_bins_abundances)
ax = gca;
ax.FontSize = 25;
% xlabel(sprintf('Unintended intron 5'' end\n(relative to intended 5'' end)'),'fontsize',18)
% ylabel('Count','fontsize',16)

subplot(1,2,2)
bar(bins,end_bins_abundances)
ax = gca;
ax.FontSize = 25;
% xlabel(sprintf('Unintended intron 3'' end\n(relative to intended 3'' end)'),'fontsize',18)
% ylabel('Count','fontsize',16)

export_fig(sprintf('%sA - cryptic_introns.png',Figures_str),'-png','-r100','-transparent');

%% Figure 5B
tmp_splicing_lib_tbl = splicing_lib_tbl([cellfind(splicing_lib_tbl.type,'endogenous'); cellfind(splicing_lib_tbl.type,'orthologous')],:);
consensus_inds=cellfind(tmp_splicing_lib_tbl.SS_5,'GTATGT');
consensus_inds = intersect(consensus_inds,cellfind(tmp_splicing_lib_tbl.branch_seq,'TACTAAC'));
tmp_tbl = tmp_splicing_lib_tbl(consensus_inds,:);

inds_AG_3 = [cellfind(tmp_tbl.SS_3,'AAG'); cellfind(tmp_tbl.SS_3,'CAG'); cellfind(tmp_tbl.SS_3,'TAG')];

inds_AAG = cellfind(tmp_tbl.SS_3,'AAG');
inds_CAG = cellfind(tmp_tbl.SS_3,'CAG');
inds_TAG = cellfind(tmp_tbl.SS_3,'TAG');
figure('units','centimeters','outerposition',[2 2 24 13.5])
% myviolinplot_mean({tmp_tbl.cryptic_splicing_eff_log(inds_AG_3),tmp_tbl.cryptic_splicing_eff_log(setdiff(1:size(tmp_tbl,1),inds_AG_3))})
% myviolinplot_mean({tmp_tbl.cryptic_splicing_eff_log(inds_AG_3),tmp_tbl.cryptic_splicing_eff_log(inds_AAG)})
myviolinplot_mean({tmp_tbl.cryptic_splicing_eff_log(inds_AAG),tmp_tbl.cryptic_splicing_eff_log(inds_CAG),tmp_tbl.cryptic_splicing_eff_log(inds_TAG),tmp_tbl.cryptic_splicing_eff_log(setdiff(1:size(tmp_tbl,1),inds_AG_3))})
% myviolinplot_mean({tmp_tbl.splicing_eff_median(inds_AAG),tmp_tbl.splicing_eff_median(inds_CAG),tmp_tbl.splicing_eff_median(inds_TAG),tmp_tbl.splicing_eff_median(setdiff(1:size(tmp_tbl,1),inds_AG_3))})
ax = gca;
ylabel('Cryptic splicing efficiency')
xlabel('Original 3'' splice site')
ax.XTickLabel = {'AAG','CAG','TAG','Other'};
ax.FontSize = 20;

export_fig(sprintf('%scryptic_intended_3SS_endogenous.png',Figures_str),'-png','-r100','-transparent');

%% Figure 5C

figure('units','centimeters','outerposition',[2 2 32 18])
myviolinplot_mean({splicing_lib_tbl.ss3_dG(inds_10), splicing_lib_tbl.ss3_dG(inds_11), splicing_lib_tbl.ss3_dG(inds_01)})
get_T_p_val_mat({splicing_lib_tbl.ss3_dG(inds_10), splicing_lib_tbl.ss3_dG(inds_11), splicing_lib_tbl.ss3_dG(inds_01)})
% ylabel('SS3 \DeltaG','fontsize',26)
ax=gca;
% ax.XTickLabel={sprintf('Only\ndesigned intron'),'Both introns','Only cryptic intron'};
ax.FontSize = 25;

export_fig(sprintf('%sB - cryptic_introns_dG.png',Figures_str),'-png','-r100','-transparent');

%% Figure 5D

a_ratio_norm=round(splicing_lib_cryptic_tbl.a_ratio*100,0);
a_ss3_mat=[];
for i=1:length(a_ratio_norm)
    a_ss3_mat = [a_ss3_mat; repmat(splicing_lib_cryptic_tbl.a_ss3_seq{i},a_ratio_norm(i),1)];
end

Alphabet='ACGT';
CT=cbrewer('qual','Set1',4);
CT_tmp=cbrewer('qual','Set2',8);
CT=[CT; CT_tmp(6,:)];
colors = cell(6,2);
    colors(1,:) = {'A', CT(3,:)};
    colors(2,:) = {'C', CT(2,:)};
    colors(3,:) = {'G', CT(5,:)};
    colors(4,:) = {'T', CT(1,:)};
    colors(5,:) = {'U', CT(1,:)};
    colors(6,:) = {'', CT(4,:)};
    
figure('units','centimeters','outerposition',[2 2 32 18])
a_ss3_mat(a_ss3_mat=='T')='U';
SeqLogoFig_fixedH(a_ss3_mat,'colors',colors,'CUTOFF',0);
% f=gcf;
% pos_logo=f.Position;
% W5_logo=pos_logo(3);
ax = gca;
ax.XTick=[];

export_fig(sprintf('%sD - ass3_logo.png',Figures_str),'-png','-r100','-transparent');

%% Figure 5E

splicing_lib_tbl.downstream_ss3_motif=nan*ones(size(splicing_lib_tbl,1),1);
for i = 1:size(splicing_lib_tbl,1)
    if isempty(splicing_lib_tbl.SS3_inds{i})
        continue
    end
    inds1=strfind(splicing_lib_tbl.seq{i}(42+splicing_lib_tbl.SS3_inds{i}(end):end),'TAG');
    inds2=strfind(splicing_lib_tbl.seq{i}(42+splicing_lib_tbl.SS3_inds{i}(end):end),'CAG');
    inds3=strfind(splicing_lib_tbl.seq{i}(42+splicing_lib_tbl.SS3_inds{i}(end):end),'AAG');
    ind = min([inds1, inds2, inds3]);
    if ~isempty(ind)
        splicing_lib_tbl.downstream_ss3_motif(i)=ind+2;
    end
end

BCs=unique(splicing_lib_cryptic_tbl.BC);
splicing_lib_cryptic_tbl.a_is_max=zeros(size(splicing_lib_cryptic_tbl,1),1);
splicing_lib_cryptic_tbl.a_is_first_ss3_motif=zeros(size(splicing_lib_cryptic_tbl,1),1);
n=1;
for i = 1:length(BCs)
    if ~mod(i,100)
        disp(i)
    end
    ind=cellfind(splicing_lib_tbl.BC,BCs{i});
    if isempty(ind)
        continue
    end
    inds=cellfind(splicing_lib_cryptic_tbl.BC,BCs{i});
    aSS3_tmp=cellfun(@(x) x(1),splicing_lib_cryptic_tbl.a_ss3_inds(inds))-splicing_lib_tbl.SS3_inds{ind}(1);
    splicing_lib_cryptic_tbl.a_ss3_shift(inds)=aSS3_tmp;
    first_ss3_inds=inds(find(abs(aSS3_tmp-splicing_lib_tbl.downstream_ss3_motif(ind))<3));
    splicing_lib_cryptic_tbl.a_is_first_ss3_motif(first_ss3_inds)=1;   
    alternative_isoforms_num(i)=length(inds);
    max_inds=find(splicing_lib_cryptic_tbl.a_ratio(inds)==max(splicing_lib_cryptic_tbl.a_ratio(inds)));
    splicing_lib_cryptic_tbl.a_is_max(inds(max_inds))=1;
    for j=1:length(max_inds)
        aSS3_seq=splicing_lib_cryptic_tbl.a_ss3_seq{inds(j)};
        if strcmp(aSS3_seq,'CAG') ||strcmp(aSS3_seq,'TAG') || strcmp(aSS3_seq,'AAG')
            is_ss3_motif(n)=1;
        else
            is_ss3_motif(n)=0;
        end
        downstream_ss3(n)=splicing_lib_tbl.downstream_ss3_motif(ind);
        aSS3(n)=aSS3_tmp(j);
        n=n+1;
    end
end

inds=cellfind(splicing_lib_cryptic_tbl.a_ss3_seq,'CAG');
inds = [inds; cellfind(splicing_lib_cryptic_tbl.a_ss3_seq,'TAG')];
inds = [inds; cellfind(splicing_lib_cryptic_tbl.a_ss3_seq,'AAG')];

motif_downstream_inds=intersect(inds,find(splicing_lib_cryptic_tbl.a_ss3_shift>0));
motif_fraction_downstream=length(motif_downstream_inds)/sum(splicing_lib_cryptic_tbl.a_ss3_shift>0)
motif_1st_fraction_downstream=sum(splicing_lib_cryptic_tbl.a_is_first_ss3_motif(motif_downstream_inds))/length(motif_downstream_inds)

a_ratio_norm=round(splicing_lib_cryptic_tbl.a_ratio*100,0);
a_ss3_mat=[];
is_first = [];
is_consensus = [];
for i=1:length(a_ratio_norm)
    a_ss3_mat = [a_ss3_mat; repmat(splicing_lib_cryptic_tbl.a_ss3_seq{i},a_ratio_norm(i),1)];
    is_first = [is_first; repmat(splicing_lib_cryptic_tbl.a_is_first_ss3_motif(i),a_ratio_norm(i),1)];
    if strcmp(splicing_lib_cryptic_tbl.a_ss3_seq{i},'TAG') | strcmp(splicing_lib_cryptic_tbl.a_ss3_seq{i},'CAG') | strcmp(splicing_lib_cryptic_tbl.a_ss3_seq{i},'AAG')
        is_consensus = [is_consensus; ones(a_ratio_norm(i),1)];
    else
        is_consensus = [is_consensus; zeros(a_ratio_norm(i),1)];
    end
end
motif_1st_fraction_downstream2=sum(is_first.*is_consensus)/sum(is_consensus)

tmp_p = splicing_lib_cryptic_tbl(motif_downstream_inds,:);
tmp_n = splicing_lib_cryptic_tbl(setdiff(1:size(splicing_lib_cryptic_tbl,1),motif_downstream_inds),:);
figure
myviolinplot_mean({(tmp_p.a_ratio(tmp_p.a_is_first_ss3_motif==1)), (tmp_n.a_ratio(tmp_n.a_is_first_ss3_motif==0))})
[~,p]=ttest2((tmp_n.a_ratio(tmp_n.a_is_first_ss3_motif==0)),(tmp_p.a_ratio(tmp_p.a_is_first_ss3_motif==1)))
% ylabel('3''SS alterntively spliced ratio','fontsize',16)
ax = gca;
% ax.XTickLabel = {'1^{st} downstream 3''SS motif','No 3''SS motif'};
ax.FontSize = 25;


export_fig(sprintf('%sE - unintended_introns_3ss_motif.png',Figures_str),'-png','-r100','-transparent');

%% Figure 5F

load('yeast_genome_introns.mat')

S=fastaread('orf_genomic_1000_all.fasta');
S_genes=struct2table(S);

organism = 'scer';

gene_name_s=cellfun(@(x) x(1: find(x==' ',1)-1),S_genes.Header,'un',0);
[~,inds1,inds2] = intersect(gene_name_s,gene_name);
S_genes=S_genes(inds1,:);
S_genes.gene_name=gene_name_s(inds1);

coordinates_folder='../Data/SI_Hooks_et_al/coordinates/';
genome_folder='../Data/complete genomes/';
files=dir(genome_folder);
files={files.name};
genome_file_ind=strmatch(organism,files);

S=fastaread([genome_folder files{genome_file_ind}]);
S_genome=struct2table(S);

C=import_hooks_coordinates([coordinates_folder organism '.gff']);
C=struct2table(C);

% load('ribosomal_genes.mat')
% [~,inds,~]=intersect(C.parent_gene,non_ribo_genes);
% C=C(inds,:);

ss_1='TAG';
ss_2='CAG';
ss_3='AAG';
% ss_1='AGT';
% ss_2='ACG';

W=300;
upstream_3ss_mat=zeros(size(C,1),W);
downstream_3ss_mat=zeros(size(C,1),W);


for i=1:size(C,1)
    
    chr_seq=upper(S_genome.Sequence{C.chr(i)});
    if C.strand{i}=='+'
        downstream_seq=chr_seq(C.end_ind(i)+1:C.end_ind(i)+W);
        upstream_seq=chr_seq(C.end_ind(i)-3-W+1:C.end_ind(i)-3);
    else
        downstream_seq=seqrcomplement(chr_seq(C.end_ind(i)-W:C.end_ind(i)-1));
        upstream_seq=seqrcomplement(chr_seq(C.end_ind(i)+3:C.end_ind(i)+2+W));
    end
    downstream_3ss_mat(i,strfind(downstream_seq,ss_2))=1;
    downstream_3ss_mat(i,strfind(downstream_seq,ss_1))=1;
    upstream_3ss_mat(i,strfind(upstream_seq,ss_2))=1;
    upstream_3ss_mat(i,strfind(upstream_seq,ss_1))=1;
    
end


M=1e3;
upstream_rnd_mean_mat=zeros(M,W);
downstream_rnd_mean_mat=zeros(M,W);
inds=find(cellfun(@length,S_genes.Sequence)>(2000+2*W+10));
for n=1:M
    upstream_3ss_mat_rnd=zeros(size(C,1),W);
    downstream_3ss_mat_rnd=zeros(size(C,1),W);
    gene_ind=inds(randi(length(inds),size(C,1),1));
    for i=1:size(C,1)
        gene_seq=upper(S_genes.Sequence{gene_ind(i)}(1001:end-1000));
        ind=randi([W+4 length(gene_seq)-W]);
        downstream_seq=gene_seq(ind+1:ind+W);
        upstream_seq=gene_seq(ind-3-W+1:ind-3);
        downstream_3ss_mat_rnd(i,strfind(downstream_seq,ss_2))=1;
        downstream_3ss_mat_rnd(i,strfind(downstream_seq,ss_1))=1;
        upstream_3ss_mat_rnd(i,strfind(upstream_seq,ss_2))=1;
        upstream_3ss_mat_rnd(i,strfind(upstream_seq,ss_1))=1;
    end
    downstream_rnd_mean_mat(n,:)=mean(downstream_3ss_mat_rnd);
    upstream_rnd_mean_mat(n,:)=mean(upstream_3ss_mat_rnd);
end


figure('units','centimeters','outerposition',[2 2 32 13.5])
% subplot(1,2,1)
hold on
% h1=plot(-W:W-2,[upstream_rnd_mean_mat(:,1:end-2)'; nan*ones(3,M) ;downstream_rnd_mean_mat(:,1:end-2)'],'-.','Color',[0,0,0]+0.8);

h3=plot(-W:W-2,[mean(upstream_3ss_mat(:,1:end-2)), nan*ones(1,3),mean(downstream_3ss_mat(:,1:end-2))],'k','linewidth',1.5);
h2=plot(-W:W-2,[mean(upstream_rnd_mean_mat(:,1:end-2)), nan*ones(1,3),mean(downstream_rnd_mean_mat(:,1:end-2))] ,'r-.','linewidth',3);
% h4=stem(-2:0,ones(1,3));
CT = cbrewer('qual','Set1',9);
h4=bar(-1,1);
h4.BarWidth=3;
h4.EdgeColor='none';
h4.FaceColor=CT(5,:);
h4.FaceAlpha=0.8;
ax = gca;
ax.FontSize = 25;
% xlabel('Positiion (relative to intron''s 3'' end)','fontsize',16)
% ylabel('Frequency of consensus 3''SS','fontsize',16)
ylim([0 0.06])
h=legend([h3,h2,h4],{'Intron containing genes','Mean of randomized set','3'' splice site'});
h.Location='best outside';
h.FontSize=16;
        
export_fig(sprintf('%sF - cerevisiae_3ss_avoidance.png',Figures_str),'-png','-r100','-transparent');


