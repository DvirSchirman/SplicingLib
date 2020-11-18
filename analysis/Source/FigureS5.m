clear
close all
clc

addpath('../Data/');
addpath(genpath('Functions/'));

Figures_str = '../Figures/FigureS5/';
if ~exist(Figures_str)
    mkdir(Figures_str)
end

coordinates_folder='../Data/SI_Hooks_et_al/coordinates/';
genome_folder='../Data/complete genomes/';
CDS_folder='../Data/CDS/';
files=dir(genome_folder);
files={files.name};
organisms={'cgla';'ncas';'ndai';'kafr';'tbla';'ecym';'zrou';'tdel';'kthe';'klac'};
species_list={'C.glabrata';'N.castellii';'N.dairenensis';'K.africana';'T.blattae';'E.cymbalariae';'Z.rouxii';'T.delbrueckii';'K.thermotolerans';'K.lactis'};

ss_1='TAG';
ss_2='CAG';
ss_3='AAG';
W=200;
M=100;
W_avoid=30;

CT = cbrewer('qual','Set1',4);
figure('units','centimeters','outerposition',[1 1 40 22.5])
for k=1:length(organisms)
    genome_file_ind=strmatch(organisms{k},files);
    S=fastaread([genome_folder files{genome_file_ind}]);
    S_genome=struct2table(S);
    
    S=fastaread([CDS_folder organisms{k} '.fasta']);
    S_genes=struct2table(S);

    C=import_hooks_coordinates([coordinates_folder organisms{k} '.gff']);
    C=struct2table(C);
    
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
        
%         downstream_3ss_mat(i,strfind(downstream_seq,ss_3))=1;
%         upstream_3ss_mat(i,strfind(upstream_seq,ss_3))=1;
    end
    
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
    %         gene_seq=upper(S_genes.Sequence{gene_ind(i)});
    %         ind=randi([1001 1200]);
            downstream_seq=gene_seq(ind+1:ind+W);
            upstream_seq=gene_seq(ind-3-W+1:ind-3);
            downstream_3ss_mat_rnd(i,strfind(downstream_seq,ss_2))=1;
            downstream_3ss_mat_rnd(i,strfind(downstream_seq,ss_1))=1;
            upstream_3ss_mat_rnd(i,strfind(upstream_seq,ss_2))=1;
            upstream_3ss_mat_rnd(i,strfind(upstream_seq,ss_1))=1;
            
%             downstream_3ss_mat_rnd(i,strfind(downstream_seq,ss_3))=1;
%             upstream_3ss_mat_rnd(i,strfind(upstream_seq,ss_3))=1;
        end
        downstream_rnd_mean_mat(n,:)=mean(downstream_3ss_mat_rnd);
        upstream_rnd_mean_mat(n,:)=mean(upstream_3ss_mat_rnd);
    end
    
    p_window_upstream(k)=sum(mean(mean(upstream_3ss_mat(:,end-1-W_avoid:end-2)))>=mean(upstream_rnd_mean_mat(:,end-1-W_avoid:end-2),2))/M;
    p_window_downstream(k)=sum(mean(mean(downstream_3ss_mat(:,1:W_avoid)))>=mean(downstream_rnd_mean_mat(:,1:W_avoid),2))/M;
    
    subplot(5,2,k)
    hold on
    plot(-W:W-2,[mean(upstream_3ss_mat(:,1:end-2)), nan*ones(1,3),mean(downstream_3ss_mat(:,1:end-2))],'k','linewidth',1.5)
    plot(-W:W-2,[mean(upstream_rnd_mean_mat(:,1:end-2)), nan*ones(1,3),mean(downstream_rnd_mean_mat(:,1:end-2))],'r-.','linewidth',1.5)
    h4=bar(-1,1);
    h4.BarWidth=3;
    h4.EdgeColor='none';
    h4.FaceAlpha=0.8;
    % plot(1/32*ones(1,size(upstream_3ss_mat,2)),'--')
    ylim([0 0.099])
    if k==5
        xlabel(['\it{' species_list{k} '}'],'fontsize',14,'Color',CT(4,:))
    elseif ismember(k,6:9)
        xlabel(['\it{' species_list{k} '}'],'fontsize',14,'Color','r')
    else
        xlabel(['\it{' species_list{k} '}'],'fontsize',14,'Color','k')
    end
end


export_fig(sprintf('%scryptic_introns.png',Figures_str),'-png','-r500','-transparent');
