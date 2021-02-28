clear
close all
clc

warning('off')

addpath('../Data/');
addpath(genpath('Functions/'));

load('SplicingLib_db.mat')
load('dG_sliding.mat')

[~,inds1,inds2]=intersect(splicing_lib_tbl.BC,dG_sliding.BC);
splicing_lib_tbl=splicing_lib_tbl(inds1,:);
dG_sliding=dG_sliding(inds2,:);

%%

get_T = @(x) (sum(x=='T'))/length(x);
get_GC =@(x) (sum(x=='C')+sum(x=='G'))./length(x);

splicing_lib_cryptic_tbl = table();

THR=0.01;
for i = 1:size(splicing_lib_tbl,1)
    if splicing_lib_tbl.cryptic_splicing_eff_log(i) < THR
        continue;
    end
    for j = 1:length(splicing_lib_tbl.cryptic_inds_abundance_log{i})
        if abs(splicing_lib_tbl.cryptic_start_inds_log{i}(j))<=3
            splicing_lib_tbl.cryptic_start_inds_log{i}(j)=0;
        end
        if abs(splicing_lib_tbl.cryptic_end_inds_log{i}(j))<=3
            splicing_lib_tbl.cryptic_end_inds_log{i}(j)=0;
        end        
    end
    inds_tuple=arrayfun(@(x,y) [num2str(x) ',' num2str(y)],splicing_lib_tbl.cryptic_start_inds_log{i},splicing_lib_tbl.cryptic_end_inds_log{i},'un',0);
    [tmp,inds]=unique(inds_tuple);
    if length(tmp)<length(inds_tuple)
        start_inds = splicing_lib_tbl.cryptic_start_inds_log{i}(inds);
        end_inds = splicing_lib_tbl.cryptic_end_inds_log{i}(inds);
        for k=1:length(inds)
            abundances(k)=sum(splicing_lib_tbl.cryptic_inds_abundance_log{i}(cellfind(inds_tuple,tmp{k})));
        end
        splicing_lib_tbl.cryptic_start_inds_log{i} = start_inds;
        splicing_lib_tbl.cryptic_end_inds_log{i} = end_inds;
        splicing_lib_tbl.cryptic_inds_abundance_log{i} = abundances;
        
        clear start_inds end_inds abundances
    end
end

n=1;
for i = 1:size(splicing_lib_tbl,1)
    if splicing_lib_tbl.cryptic_splicing_eff_log(i) < THR
        continue;
    end
    for j = 1:length(splicing_lib_tbl.cryptic_inds_abundance_log{i})
        if splicing_lib_tbl.cryptic_inds_abundance_log{i}(j) > 4
            
            
            splicing_lib_cryptic_tbl.BC{n}=splicing_lib_tbl.BC{i};
            splicing_lib_cryptic_tbl.type{n}=splicing_lib_tbl.type{i};
            splicing_lib_cryptic_tbl.intended_splicing{n}=splicing_lib_tbl.splicing_eff_median(i);
            splicing_lib_cryptic_tbl.a_ratio(n)=splicing_lib_tbl.cryptic_inds_abundance_log{i}(j)./splicing_lib_tbl.total_abundance_log(i);
            splicing_lib_cryptic_tbl.a_ss5_inds{n}=splicing_lib_tbl.SS5_inds{i}+splicing_lib_tbl.cryptic_start_inds_log{i}(j);
            splicing_lib_cryptic_tbl.a_ss5_seq{n}=splicing_lib_tbl.seq{i}(42+splicing_lib_tbl.SS5_inds{i}+splicing_lib_tbl.cryptic_start_inds_log{i}(j));
            splicing_lib_cryptic_tbl.a_ss3_inds{n}=splicing_lib_tbl.SS3_inds{i}+splicing_lib_tbl.cryptic_end_inds_log{i}(j)-2;
            splicing_lib_cryptic_tbl.a_ss3_seq{n}=splicing_lib_tbl.seq{i}(42+splicing_lib_tbl.SS3_inds{i}+splicing_lib_tbl.cryptic_end_inds_log{i}(j)-2);
            splicing_lib_cryptic_tbl.a_branch_seq{n}=splicing_lib_tbl.branch_seq{i};
            splicing_lib_cryptic_tbl.a_intron_seq{n}=splicing_lib_tbl.seq{i}(42+splicing_lib_tbl.SS5_inds{i}(1)+splicing_lib_tbl.cryptic_start_inds_log{i}(j):...
                                                            (42+splicing_lib_tbl.SS3_inds{i}(end)+splicing_lib_tbl.cryptic_end_inds_log{i}(j)-2));
            splicing_lib_cryptic_tbl.a_inron_len(n)=length(splicing_lib_cryptic_tbl.a_intron_seq{n});
            if ~isempty(splicing_lib_tbl.branch_inds{i})
                splicing_lib_cryptic_tbl.a_branch_to_3(n)=42+splicing_lib_tbl.SS3_inds{i}(1)+splicing_lib_tbl.cryptic_end_inds_log{i}(j)-2-splicing_lib_tbl.branch_inds{i}(end);
            end
            splicing_lib_cryptic_tbl.a_GC(n)=get_GC(splicing_lib_cryptic_tbl.a_intron_seq{n});
            splicing_lib_cryptic_tbl.a_polyT(n)=get_T(splicing_lib_cryptic_tbl.a_intron_seq{n}(end-22:end-3));
            if 27+splicing_lib_tbl.SS3_inds{i}(2)+splicing_lib_tbl.cryptic_end_inds_log{i}(j)<=200
                splicing_lib_cryptic_tbl.a_ss3_dG(n)=dG_sliding.dG{i}(27+splicing_lib_tbl.SS3_inds{i}(2)+splicing_lib_tbl.cryptic_end_inds_log{i}(j)-2);
            end
            
            n=n+1;
        end
    end
end

%%
ctrl_inds = find(cell2mat(cellfun(@(x) ~isempty(strfind(x,'control')),splicing_lib_cryptic_tbl.type,'un',0)));
% ctrl_inds(find(cell2mat(cellfun(@(x) ~isempty(strfind(x,'intron')),splicing_lib_tbl.type(ctrl_inds),'un',0))))=[];
% ctrl_inds = [ctrl_inds; cellfind(splicing_lib_tbl.orth_organism,'agos')];

ctrl_splicing_lib_cryptic_tbl=splicing_lib_cryptic_tbl(ctrl_inds,:);
splicing_lib_cryptic_tbl=splicing_lib_cryptic_tbl(setdiff(1:size(splicing_lib_cryptic_tbl,1),ctrl_inds),:);

inds = find(cell2mat(cellfun(@(x) ~isempty(strfind(x,'intron')),ctrl_splicing_lib_cryptic_tbl.type,'un',0)));
ctrl_splicing_lib_cryptic_tbl(inds,:)=[];

%%
inds = find(cellfun(@isempty, splicing_lib_tbl.SS5_inds));    
splicing_lib_tbl(inds,:)=[];
inds = find(cellfun(@isempty, splicing_lib_tbl.SS3_inds));    
splicing_lib_tbl(inds,:)=[];

ctrl_inds = find(cell2mat(cellfun(@(x) ~isempty(strfind(x,'control')),splicing_lib_tbl.type,'un',0)));
% ctrl_inds(find(cell2mat(cellfun(@(x) ~isempty(strfind(x,'intron')),splicing_lib_tbl.type(ctrl_inds),'un',0))))=[];
ctrl_inds = [ctrl_inds; cellfind(splicing_lib_tbl.orth_organism,'agos')];
ctrl_tbl=splicing_lib_tbl(ctrl_inds,:);

splicing_lib_tbl=splicing_lib_tbl(setdiff(1:size(splicing_lib_tbl,1),ctrl_inds),:);

inds = find(cell2mat(cellfun(@(x) ~isempty(strfind(x,'intron')),ctrl_tbl.type,'un',0)));
ctrl_tbl(inds,:)=[];

splicing_lib_cryptic_tbl = splicing_lib_cryptic_tbl(isfinite(splicing_lib_cryptic_tbl.a_ratio),:);
ctrl_splicing_lib_cryptic_tbl = ctrl_splicing_lib_cryptic_tbl(isfinite(ctrl_splicing_lib_cryptic_tbl.a_ratio),:);

save('../Data/unintended_splicing_tbl.mat', 'splicing_lib_cryptic_tbl','ctrl_splicing_lib_cryptic_tbl')

            