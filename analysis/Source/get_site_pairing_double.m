clear
close all
clc

addpath('../Data/');
addpath(genpath('Functions/'));

load('SplicingLib_db.mat')
load('lib_structures_tbl.mat')
load('double_sp_eff_tbl.mat')

[~,inds1,inds2]=intersect(splicing_lib_tbl.BC,double_sp_eff_tbl.BC);
double_sp_eff_tbl=double_sp_eff_tbl(inds2,:);
splicing_lib_tbl=splicing_lib_tbl(inds1,:);
transcript_structure_tbl=transcript_structure_tbl(inds1,:);


%%

for i=1:size(double_sp_eff_tbl,1)
    ss5_inds=double_sp_eff_tbl.SS5_1_inds{i}+31;
    ss5_vec=ones(1,6);
    ss5_struct=transcript_structure_tbl.Structure{i}(ss5_inds:ss5_inds+5);
    ss5_vec(ss5_struct=='.')=0;
    double_sp_eff_tbl.ss5_1_pairing{i}=ss5_vec;
    
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
        double_sp_eff_tbl.ss5_1_stem_len(i)=open_first3-open_first5-1;
        if sum(ss5_struct=='(')==6
            double_sp_eff_tbl.ss5_1_pairing_dir(i)=1;
        elseif sum(ss5_struct==')')==6
            double_sp_eff_tbl.ss5_1_pairing_dir(i)=-1;
        else
            double_sp_eff_tbl.ss5_1_pairing_dir(i)=0;
        end        
    else
        double_sp_eff_tbl.ss5_1_stem_len(i)=0;
        double_sp_eff_tbl.ss5_1_pairing_dir(i)=0;
    end
    
%     branch_inds=double_sp_eff_tbl.branch_ind_1(i)+31:double_sp_eff_tbl.branch_ind_1(i)+31+6;
    branch_inds=double_sp_eff_tbl.branch_1_inds{i}+31;
    branch_vec=ones(1,7);
    branch_struct=transcript_structure_tbl.Structure{i}(branch_inds:branch_inds+5);
    branch_vec(branch_struct=='.')=0;
    double_sp_eff_tbl.branch_1_pairing{i}=branch_vec;
    
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
        double_sp_eff_tbl.branch_1_stem_len(i)=open_first3-open_first5-1;
        if sum(branch_struct=='(')==7
            double_sp_eff_tbl.branch_1_pairing_dir(i)=1;
        elseif sum(branch_struct==')')==7
            double_sp_eff_tbl.branch_1_pairing_dir(i)=-1;
        else
            double_sp_eff_tbl.branch_1_pairing_dir(i)=0;
        end        
    else
        double_sp_eff_tbl.branch_1_stem_len(i)=0;
        double_sp_eff_tbl.branch_1_pairing_dir(i)=0;
    end
    
    ss3_inds=double_sp_eff_tbl.SS3_1_inds{i}+31;
    ss3_vec=ones(1,3);
    ss3_struct=transcript_structure_tbl.Structure{i}(ss3_inds:ss3_inds+2);
    ss3_vec(ss3_struct=='.')=0;
    double_sp_eff_tbl.ss3_1_pairing{i}=ss3_vec;
    
    if sum(ss3_vec)==3
        extended_struct=transcript_structure_tbl.Structure{i}(ss3_inds-30:ss3_inds+2+30);
        open_inds=find(extended_struct=='.');
        tmp=31-open_inds;
        tmp(tmp<0)=nan;
        [~,open_first5]=min(tmp);
        open_first5=open_inds(open_first5);
        tmp=open_inds-33;
        tmp(tmp<0)=nan;
        [~,open_first3]=min(tmp);
        open_first3=open_inds(open_first3);
        double_sp_eff_tbl.ss3_1_stem_len(i)=open_first3-open_first5-1;
        if sum(ss3_struct=='(')==3
            double_sp_eff_tbl.ss3_1_pairing_dir(i)=1;
        elseif sum(ss3_struct==')')==3
            double_sp_eff_tbl.ss3_1_pairing_dir(i)=-1;
        else
            double_sp_eff_tbl.ss3_1_pairing_dir(i)=0;
        end        
    else
        double_sp_eff_tbl.ss3_1_stem_len(i)=0;
        double_sp_eff_tbl.ss3_1_pairing_dir(i)=0;
    end
    
    
    ss5_inds=double_sp_eff_tbl.SS5_2_inds{i}+31;
    ss5_vec=ones(1,6);
    ss5_struct=transcript_structure_tbl.Structure{i}(ss5_inds:ss5_inds+5);
    ss5_vec(ss5_struct=='.')=0;
    double_sp_eff_tbl.ss5_2_pairing{i}=ss5_vec;
    
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
        double_sp_eff_tbl.ss5_2_stem_len(i)=open_first3-open_first5-1;
        if sum(ss5_struct=='(')==6
            double_sp_eff_tbl.ss5_2_pairing_dir(i)=1;
        elseif sum(ss5_struct==')')==6
            double_sp_eff_tbl.ss5_2_pairing_dir(i)=-1;
        else
            double_sp_eff_tbl.ss5_2_pairing_dir(i)=0;
        end        
    else
        double_sp_eff_tbl.ss5_2_stem_len(i)=0;
        double_sp_eff_tbl.ss5_2_pairing_dir(i)=0;
    end
    
%     branch_inds=double_sp_eff_tbl.branch_ind_1(i)+31:double_sp_eff_tbl.branch_ind_1(i)+31+6;
    branch_inds=double_sp_eff_tbl.branch_2_inds{i}+31;
    branch_vec=ones(1,7);
    branch_struct=transcript_structure_tbl.Structure{i}(branch_inds:branch_inds+5);
    branch_vec(branch_struct=='.')=0;
    double_sp_eff_tbl.branch_2_pairing{i}=branch_vec;
    
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
        double_sp_eff_tbl.branch_2_stem_len(i)=open_first3-open_first5-1;
        if sum(branch_struct=='(')==7
            double_sp_eff_tbl.branch_2_pairing_dir(i)=1;
        elseif sum(branch_struct==')')==7
            double_sp_eff_tbl.branch_2_pairing_dir(i)=-1;
        else
            double_sp_eff_tbl.branch_2_pairing_dir(i)=0;
        end        
    else
        double_sp_eff_tbl.branch_2_stem_len(i)=0;
        double_sp_eff_tbl.branch_2_pairing_dir(i)=0;
    end
    
    ss3_inds=double_sp_eff_tbl.SS3_2_inds{i}+31;
    ss3_vec=ones(1,3);
    ss3_struct=transcript_structure_tbl.Structure{i}(ss3_inds:ss3_inds+2);
    ss3_vec(ss3_struct=='.')=0;
    double_sp_eff_tbl.ss3_2_pairing{i}=ss3_vec;
    
    if sum(ss3_vec)==3
        end_inds=min([ss3_inds(end)+2+30,210]);
        extended_struct=transcript_structure_tbl.Structure{i}(ss3_inds-30:end_inds);
        open_inds=find(extended_struct=='.');
        tmp=31-open_inds;
        tmp(tmp<0)=nan;
        [~,open_first5]=min(tmp);
        open_first5=open_inds(open_first5);
        tmp=open_inds-33;
        tmp(tmp<0)=nan;
        [~,open_first3]=min(tmp);
        open_first3=open_inds(open_first3);
        double_sp_eff_tbl.ss3_2_stem_len(i)=open_first3-open_first5-1;
        if sum(ss3_struct=='(')==3
            double_sp_eff_tbl.ss3_2_pairing_dir(i)=1;
        elseif sum(ss3_struct==')')==3
            double_sp_eff_tbl.ss3_2_pairing_dir(i)=-1;
        else
            double_sp_eff_tbl.ss3_2_pairing_dir(i)=0;
        end        
    else
        double_sp_eff_tbl.ss3_2_stem_len(i)=0;
        double_sp_eff_tbl.ss3_2_pairing_dir(i)=0;
    end
end  

save('../Data/double_site_structure_pairing.mat','double_sp_eff_tbl')
% %%    
% [~,inds]=sort(splicing_lib_tbl.splicing_eff_median,'descend','MissingPlacement','last');
% figure
% subplot(1,3,1)
% imagesc(cell2mat(site_pairing_tbl.ss5_pairing(inds)))
% subplot(1,3,2)
% imagesc(cell2mat(site_pairing_tbl.branch_pairing(inds)))
% subplot(1,3,3)
% imagesc(cell2mat(site_pairing_tbl.ss3_pairing(inds)))
