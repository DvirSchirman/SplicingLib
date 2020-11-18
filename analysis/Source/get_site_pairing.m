clear
close all
clc

addpath('../Data/');
addpath(genpath('Functions/'));

load('SplicingLib_db.mat')
load('lib_structures_tbl.mat')

[~,inds1,inds2] = intersect(splicing_lib_tbl.BC,transcript_structure_tbl.BC);
splicing_lib_tbl = splicing_lib_tbl(inds1,:);
transcript_structure_tbl = transcript_structure_tbl(inds2,:);

%%
load('double_intron_db.mat')
% inds=find(cell2mat(cellfun(@(x) ~isempty(strfind(x,'joined')),double_tbl.type,'un',0)));
% double_tbl=double_tbl(inds,:);
% double_tbl(cellfind(double_tbl.type,'double endogenous - joined ss control'),:)=[];
% double_tbl(:,cellfind(double_tbl.Properties.VariableNames,'SS5_inds'))=[];
% double_tbl(:,cellfind(double_tbl.Properties.VariableNames,'SS3_inds'))=[];
% 
% double_tbl.Properties.VariableNames{cellfind(double_tbl.Properties.VariableNames,'SS5_1_inds')}='SS5_inds';
% double_tbl.Properties.VariableNames{cellfind(double_tbl.Properties.VariableNames,'SS3_1_inds')}='SS3_inds';
% double_tbl.SS_5=cellfun(@(x,y) x(y+12),double_tbl.seq,double_tbl.SS5_inds,'un',0);
% double_tbl.SS_3=cellfun(@(x,y) x(y+12),double_tbl.seq,double_tbl.SS3_inds,'un',0);
% double_tbl.intron_len=cellfun(@(x,y) y(end)-x(1),double_tbl.SS5_inds,double_tbl.SS3_inds);
% 
% branch_inds=cellfun(@(x) strfind(x,'TAAC')-3-12,double_tbl.seq,'un',0);
% double_tbl(find(cellfun(@isempty,branch_inds)),:)=[];
% branch_inds(find(cellfun(@isempty,branch_inds)))=[];
% branch_from_3=cellfun(@(x,y) y(end)-x,branch_inds,double_tbl.SS3_inds,'un',0);
% for i=1:length(branch_inds)
%     if length(branch_inds{i})>1
%         [a,b]=min(abs(branch_from_3{i}-15));
%         branch_inds{i}=branch_inds{i}(b);
%         branch_from_3{i}=branch_from_3{i}(b);
%     end
% end
% double_tbl(find(cellfun(@(x) x<1, branch_inds)),:)=[];
% branch_inds(find(cellfun(@(x) x<1, branch_inds)))=[];
% 
% branch_seqs=cellfun(@(x,y) x(y+12:y+12+6),double_tbl.seq,branch_inds,'un',0);
% 
% double_tbl.branch_inds=branch_inds;
% double_tbl.branch_seq=branch_seqs;
% %%
% 
% [~,inds1,inds2]=intersect(splicing_lib_tbl.BC,double_tbl.BC);
% 
% splicing_lib_tbl.SS_5(inds1)=double_tbl.SS_5(inds2);
% splicing_lib_tbl.SS5_inds(inds1)=double_tbl.SS5_inds(inds2);
% splicing_lib_tbl.SS_3(inds1)=double_tbl.SS_3(inds2);
% splicing_lib_tbl.SS3_inds(inds1)=double_tbl.SS3_inds(inds2);
% splicing_lib_tbl.branch_inds(inds1)=double_tbl.branch_inds(inds2);
% splicing_lib_tbl.branch_seq(inds1)=double_tbl.branch_seq(inds2);
% splicing_lib_tbl.intron_len(inds1)=mat2cell(double_tbl.intron_len(inds2),ones(size(double_tbl,1),1));

%%

inds = find(cellfun(@isempty, splicing_lib_tbl.SS5_inds));    
splicing_lib_tbl(inds,:)=[];
transcript_structure_tbl(inds,:)=[];

% [~,inds1,inds2]=intersect(splicing_lib_tbl.BC,transcript_structure_tbl.BC);
% sum(cell2mat(cellfun(@(x,y) strcmp(x,y),splicing_lib_tbl.BC,transcript_structure_tbl.BC,'un',0)))
site_pairing_tbl=splicing_lib_tbl(:,cellfind(splicing_lib_tbl.Properties.VariableNames,'BC'));

for i=1:size(transcript_structure_tbl,1)
    ss5_inds=splicing_lib_tbl.SS5_inds{i}(1)+31;
    ss3_inds=splicing_lib_tbl.SS3_inds{i}(1)+31;
    branch_inds=strfind(transcript_structure_tbl.Sequence{i},splicing_lib_tbl.branch_seq{i});
%     branch_inds=splicing_lib_tbl.branch_inds{i}(1)+31;
    if isempty(branch_inds) & ~isempty(splicing_lib_tbl.branch_seq{i})
        branch_inds=strfind(transcript_structure_tbl.Sequence{i},splicing_lib_tbl.branch_seq{i}(2:end))-1;
    end
    if length(branch_inds)>1
        tmp=ss3_inds-branch_inds;
        tmp(tmp<0)=nan;
        [~,b]=min(tmp);
        branch_inds=branch_inds(b);
    end
    ss5_vec=ones(1,6);
    ss5_struct=transcript_structure_tbl.Structure{i}(ss5_inds:ss5_inds+5);
    ss5_vec(ss5_struct=='.')=0;
    site_pairing_tbl.ss5_pairing{i}=ss5_vec;
    if sum(ss5_vec)==6
        extended_struct=transcript_structure_tbl.Structure{i}(ss5_inds-30:ss5_inds+5+30);
        open_inds=find(extended_struct=='.');
        tmp=31-open_inds;
        tmp(tmp<0)=nan;
        [~,open_first5]=min(tmp);
        open_first5=open_inds(open_first5);
        tmp=open_inds-36;
        tmp(tmp<0)=nan;
        [~,open_first3]=min(tmp);
        open_first3=open_inds(open_first3);
        site_pairing_tbl.ss5_stem_len(i)=open_first3-open_first5-1;
        if sum(ss5_struct=='(')==6
            site_pairing_tbl.ss5_pairing_dir(i)=1;
        elseif sum(ss5_struct==')')==6
            site_pairing_tbl.ss5_pairing_dir(i)=-1;
        else
            site_pairing_tbl.ss5_pairing_dir(i)=0;
        end        
    else
        site_pairing_tbl.ss5_stem_len(i)=0;
        site_pairing_tbl.ss5_pairing_dir(i)=0;
    end
    
    if isempty(branch_inds)
        branch_vec=nan*ones(1,7);
    else
        branch_vec=ones(1,7);
        branch_struct=transcript_structure_tbl.Structure{i}(branch_inds:branch_inds+6);
        branch_vec(branch_struct=='.')=0;
    end
    site_pairing_tbl.branch_pairing{i}=branch_vec;
    
    if isempty(branch_inds)
        site_pairing_tbl.branch_stem_len(i)=nan;
        site_pairing_tbl.branch_pairing_dir(i)=nan;
    else
        if sum(branch_vec)==7
            extended_struct=transcript_structure_tbl.Structure{i}(branch_inds-30:branch_inds+6+30);
            open_inds=find(extended_struct=='.');
            tmp=31-open_inds;
            tmp(tmp<0)=nan;
            [~,open_first5]=min(tmp);
            open_first5=open_inds(open_first5);
            tmp=open_inds-37;
            tmp(tmp<0)=nan;
            [~,open_first3]=min(tmp);
            open_first3=open_inds(open_first3);
            site_pairing_tbl.branch_stem_len(i)=open_first3-open_first5-1;
            if sum(branch_struct=='(')==7
                site_pairing_tbl.branch_pairing_dir(i)=1;
            elseif sum(branch_struct==')')==7
                site_pairing_tbl.branch_pairing_dir(i)=-1;
            else
                site_pairing_tbl.branch_pairing_dir(i)=0;
            end        
        else
            site_pairing_tbl.branch_stem_len(i)=0;
            site_pairing_tbl.branch_pairing_dir(i)=0;
        end
    end
    
    ss3_vec=ones(1,3);
    ss3_struct=transcript_structure_tbl.Structure{i}(ss3_inds:ss3_inds+2);
    ss3_vec(ss3_struct=='.')=0;
    site_pairing_tbl.ss3_pairing{i}=ss3_vec;
    
    if sum(ss3_vec)==3
        extended_struct=transcript_structure_tbl.Structure{i}(ss3_inds-30:ss3_inds+2+25);
        open_inds=find(extended_struct=='.');
        tmp=31-open_inds;
        tmp(tmp<0)=nan;
        [~,open_first5]=min(tmp);
        open_first5=open_inds(open_first5);
        tmp=open_inds-33;
        tmp(tmp<0)=nan;
        [~,open_first3]=min(tmp);
        open_first3=open_inds(open_first3);
        site_pairing_tbl.ss3_stem_len(i)=open_first3-open_first5-1;
        if sum(ss3_struct=='(')==3
            site_pairing_tbl.ss3_pairing_dir(i)=1;
        elseif sum(ss3_struct==')')==3
            site_pairing_tbl.ss3_pairing_dir(i)=-1;
        else
            site_pairing_tbl.ss3_pairing_dir(i)=0;
        end        
    else
        site_pairing_tbl.ss3_stem_len(i)=0;
        site_pairing_tbl.ss3_pairing_dir(i)=0;
    end
    
end
    
save('../Data/site_structure_pairing.mat','site_pairing_tbl')
