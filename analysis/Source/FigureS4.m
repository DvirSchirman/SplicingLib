clear
close all
clc

addpath('../Data/');
addpath(genpath('Functions/'));

Figures_str = '../Figures/FigureS4/';
if ~exist(Figures_str)
    mkdir(Figures_str)
end

get_GC = @(x) (sum(x=='G')+sum(x=='C'))/length(x);
get_T = @(x) (sum(x=='T'))/length(x);

load('SplicingLib_db.mat')
load('orthologous_introns_long.mat')

orth_tbl = struct2table(orth_intron);
orth_tbl(cellfind(orth_tbl.organism,'agos'),:)=[];
orth_tbl.len = orth_tbl.len-10;
orth_tbl.seq = cellfun(@(x) x(6:end-5),orth_tbl.seq,'un',0);
orth_tbl.ss5 = cellfun(@(x) x(1:6),orth_tbl.seq,'un',0);
tmp_tbl = splicing_lib_tbl(cellfind(splicing_lib_tbl.type,'endogenous'),:);
tmp1=tmp_tbl(tmp_tbl.BC_num==1,:);
tmp2=tmp_tbl(tmp_tbl.BC_num==4,:);
tmp2=tmp2(1:4:end,:);
tmp_tbl = [tmp1; tmp2];
tmp_tbl.intron_seq=cellfun(@(x,y,z) x(y(1)+42:z(end)+42),tmp_tbl.seq,tmp_tbl.SS5_inds,tmp_tbl.SS3_inds,'un',0);

tmp_tbl.intron_len = cell2mat(tmp_tbl.intron_len);
endogenous_tbl = table();
endogenous_tbl = orth_tbl(1,:);
endogenous_tbl(1,:)=[];
endogenous_tbl.len(1:size(tmp_tbl,1)) = tmp_tbl.intron_len;
endogenous_tbl.ss5 = tmp_tbl.SS_5;
endogenous_tbl.seq = tmp_tbl.intron_seq;
endogenous_tbl.organism = repmat({'scer'},size(endogenous_tbl,1),1);
endogenous_tbl.cervisiae_gene = tmp_tbl.intron_cerevisiae_gene;
orth_tbl = [orth_tbl; endogenous_tbl];
orth_tbl(find(orth_tbl.len>200),:)=[];
orth_tbl.GC=cellfun(@(x) get_GC(x), orth_tbl.seq);
orth_tbl.ss3=cellfun(@(x) x(end-2:end),orth_tbl.seq,'un',0);

for i=1:size(orth_tbl,1)
    branch_ind=strfind(orth_tbl.seq{i},'CTAAC')-2;
    if isempty(branch_ind)
        branch_ind=strfind(orth_tbl.seq{i},'CTAAT')-2;
        if isempty(branch_ind)
            branch_ind=strfind(orth_tbl.seq{i},'TTAAC')-2;
        end
    end
    if length(branch_ind)>1
        while orth_tbl.len(i)-branch_ind(end) < 15            
            branch_ind(end)=[];
            if isempty(branch_ind)
                break;
            end
        end
        if ~isempty(branch_ind)
            branch_ind=branch_ind(end);
        end
    end 
    branch_ind(branch_ind<10)=[];
    if ~isempty(branch_ind) 
        orth_tbl.branch{i} = orth_tbl.seq{i}(branch_ind:branch_ind+6);
        orth_tbl.branch_ind(i)=branch_ind;
    end
end
orth_tbl(find(cellfun(@isempty,orth_tbl.branch)),:)=[];
orth_tbl.BS_to_3ss = arrayfun(@(x,y) y-x-9,orth_tbl.branch_ind,orth_tbl.len);

species_list={'scer';'cgla';'ncas';'ndai';'kafr';'tbla';'ecym';'zrou';'tdel';'kthe';'klac'};
species_list_full={'S.cerevisiae';'C.glabrata';'N.castellii';'N.dairenensis';'K.africana';'T.blattae';'E.cymbalariae';'Z.rouxii';'T.delbrueckii';'K.thermotolerans';'K.lactis'};

%% Figure S4A

for i=1:length(species_list)
    inds = cellfind(orth_tbl.organism,species_list{i});
    intron_cell{i}=orth_tbl.BS_to_3ss(inds);
end
figure('units','centimeters','outerposition',[2 2 27 13.5])
myviolinplot(intron_cell)

ax = gca;
ax.XTickLabel = species_list_full;
ax.XTickLabelRotation=-45;
ax.FontSize = 25;
ylabel('BS-to-3''SS length','FontSize',20)

export_fig(sprintf('%sA- orth_BS_to_3ss.png',Figures_str),'-png','-r100','-transparent'); 

%% Figure S4B

for i=1:length(species_list)
    inds = cellfind(orth_tbl.organism,species_list{i});
    intron_cell{i}=orth_tbl.len(inds);
end
figure('units','centimeters','outerposition',[2 2 27 13.5])
myviolinplot(intron_cell)

ax = gca;
ax.XTickLabel = species_list_full;
ax.XTickLabelRotation=-45;
ax.FontSize = 25;
ylabel('Intron length','FontSize',20)

export_fig(sprintf('%sB - orth_intron_len.png',Figures_str),'-png','-r100','-transparent'); 


%% Figure S4C

for i=1:length(species_list)
    inds = cellfind(orth_tbl.organism,species_list{i});
    intron_cell{i}=orth_tbl.branch_ind(inds);
end
figure('units','centimeters','outerposition',[2 2 27 13.5])
myviolinplot(intron_cell)

ax = gca;
ax.XTickLabel = species_list_full;
ax.XTickLabelRotation=-45;
ax.FontSize = 25;
ylabel('5''SS-to-BS length','FontSize',20)
    
export_fig(sprintf('%sC - orth_5ss_to_BS.png',Figures_str),'-png','-r100','-transparent'); 

%% Figure S4D

for i=1:length(species_list)
    inds = cellfind(orth_tbl.organism,species_list{i});
    intron_cell{i}=orth_tbl.GC(inds);
end
figure('units','centimeters','outerposition',[2 2 27 13.5])
myviolinplot(intron_cell)
ylabel('%GC')
ax = gca;
ax.XTickLabel = species_list_full;
ax.XTickLabelRotation=-45;
ax.FontSize = 25;
    
export_fig(sprintf('%sD - orth_GC.png',Figures_str),'-png','-r100','-transparent'); 

%% Figure S4E

orth_tbl.polyU = cellfun(@(x) get_T(x(end-22:end-3)),orth_tbl.seq);
for i=1:length(species_list)
    inds = cellfind(orth_tbl.organism,species_list{i});
    intron_cell{i}=orth_tbl.polyU(inds);
end
figure('units','centimeters','outerposition',[2 2 27 13.5])
myviolinplot(intron_cell)

ax = gca;
ax.XTickLabel = species_list_full;
ax.XTickLabelRotation=-45;
ax.FontSize = 25;
ylabel('U - enrichment','FontSize',20)

export_fig(sprintf('%sE - orth_polyU.png',Figures_str),'-png','-r100','-transparent'); 






