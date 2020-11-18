close all
clear
clc

files = dir('alignments\CDS\*.stk');
files = {files.name}';

U2AF1_plus = {'zrou';'tdel';'agos';'ecym';'sklu';'kthe';'kwal'};
U2AF1_minus = {'scer';'smik';'skud';'suva';'cgla';'kafr';'knag';'ncas';'ndai';'klac'};
U2AF1_mut = {'tbla';'tpha';'kpol'};

n=1;
for i = 1:length(files)
    tmp=readtable(['alignments\CDS\' files{i}],'FileType','text','ReadVariableNames',0,'HeaderLines',3,'Delimiter',' ','MultipleDelimsAsOne',1);
    site_str=tmp{end-1,1};
    alignment_tbl=tmp(1:end-3,:);
    alignment_tbl.Properties.VariableNames={'id','Sequence'};
    alignment_tbl.species = cellfun(@(x) x(1:4),alignment_tbl.id,'un',0);
    lengths = cellfun(@length,alignment_tbl.Sequence);
    ss3_inds=cell2mat(strfind(site_str,'[3]'));
    ss5_inds=cell2mat(strfind(site_str,'[=5==]'));
    for j=1:length(ss3_inds)
        alignment_tbl_tmp = alignment_tbl;
        ss3_ind=ss3_inds(j);
        ss5_ind=ss5_inds(j);
        inds = find(lengths<ss3_ind);
        alignment_tbl_tmp(inds,:)=[];
        is_ss3=find(cell2mat(cellfun(@(x) x(ss3_ind+1)=='A' & x(ss3_ind+2)=='G',alignment_tbl_tmp.Sequence,'un',0)));
        alignment_tbl_tmp = alignment_tbl_tmp(is_ss3,:);
        tmp = cellfun(@(x) x(ss5_ind:ss3_ind+2),alignment_tbl_tmp.Sequence,'un',0);
        tmp = cellfun(@(x) x(x~='.' & x~='-'),tmp,'un',0);
        alignment_tbl_tmp.intron_seq=tmp;
        
        [~,plus_inds,~] = intersect(alignment_tbl_tmp.species,U2AF1_plus);
        plus_u = alignment_tbl_tmp.intron_seq(plus_inds);
        plus_u = cellfun(@(x) fliplr(x(1:end-3)),plus_u,'un',0);
        plus_entropy_u{n,1} = fliplr(get_entropy(plus_u));
        if length(ss3_inds)>1 & j==1
            plus_d = cellfun(@(x) x(ss3_ind+3:ss5_inds(2)-1),alignment_tbl.Sequence(plus_inds),'un',0);
        else
            plus_d = cellfun(@(x) x(ss3_ind+3:end),alignment_tbl.Sequence(plus_inds),'un',0);
        end
        plus_d = cellfun(@(x) x(x~='.' & x~='-'),plus_d,'un',0);
        if j==1
            if mod(ss5_ind-1,3)==0
                plus_entropy_d{n,1} = get_entropy(plus_d);
            elseif mod(ss5_ind-1,3)==1
                plus_entropy_d{n,1} = [nan, get_entropy(plus_d)];
            elseif mod(ss5_ind-1,3)==2
                plus_entropy_d{n,1} = [nan, nan, get_entropy(plus_d)];
            end
        else
            tmp=ss5_inds(1)-1+(ss5_ind-(ss3_inds(1)+3));
            if mod(tmp,3)==0
                plus_entropy_d{n,1} = get_entropy(plus_d);
            elseif mod(tmp,3)==1
                plus_entropy_d{n,1} = [nan, get_entropy(plus_d)];
            elseif mod(tmp,3)==2
                plus_entropy_d{n,1} = [nan, nan, get_entropy(plus_d)];
            end
        end
            
        [~,minus_inds,~] = intersect(alignment_tbl_tmp.species,U2AF1_minus);
        [~,mut_inds,~] = intersect(alignment_tbl_tmp.species,U2AF1_mut);
        minus_inds = [minus_inds; mut_inds];

        minus_u = alignment_tbl_tmp.intron_seq(minus_inds);
        minus_u = cellfun(@(x) fliplr(x(1:end-3)),minus_u,'un',0);
        minus_entropy_u{n,1} = fliplr(get_entropy(minus_u));
        
        if length(ss3_inds)>1 & j==1
            minus_d = cellfun(@(x) x(ss3_ind+3:ss5_inds(2)-1),alignment_tbl.Sequence(minus_inds),'un',0);
        else
            minus_d = cellfun(@(x) x(ss3_ind+3:end),alignment_tbl.Sequence(minus_inds),'un',0);
        end
        minus_d = cellfun(@(x) x(x~='.' & x~='-'),minus_d,'un',0);
        if j==1
            if mod(ss5_ind-1,3)==0
                minus_entropy_d{n,1} = get_entropy(minus_d);
            elseif mod(ss5_ind-1,3)==1
                minus_entropy_d{n,1} = [nan, get_entropy(minus_d)];
            elseif mod(ss5_ind-1,3)==2
                minus_entropy_d{n,1} = [nan, nan, get_entropy(minus_d)];
            end
        else
            tmp=ss5_inds(1)-1+(ss5_ind-(ss3_inds(1)+3));
            if mod(tmp,3)==0
                minus_entropy_d{n,1} = get_entropy(minus_d);
            elseif mod(tmp,3)==1
                minus_entropy_d{n,1} = [nan, get_entropy(minus_d)];
            elseif mod(tmp,3)==2
                minus_entropy_d{n,1} = [nan, nan, get_entropy(minus_d)];
            end
        end
        n=n+1;
    end
    
end

%%
plus_entropy_u=plus_entropy_u(cellfun(@(x) ~isempty(x),plus_entropy_u));
max_len=max(cellfun(@length,plus_entropy_u));
plus_entropy_u_mat=nan*ones(size(plus_entropy_u,1),max_len);
for i=1:size(plus_entropy_u,1)
    len=length(plus_entropy_u{i});
    plus_entropy_u_mat(i,end-len+1:end)=plus_entropy_u{i};
end

minus_entropy_u=minus_entropy_u(cellfun(@(x) ~isempty(x),minus_entropy_u));
max_len=max(cellfun(@length,minus_entropy_u));
minus_entropy_u_mat=nan*ones(size(minus_entropy_u,1),max_len);
for i=1:size(minus_entropy_u,1)
    len=length(minus_entropy_u{i});
    minus_entropy_u_mat(i,end-len+1:end)=minus_entropy_u{i};
end

plus_entropy_d=plus_entropy_d(cellfun(@(x) ~isempty(x),plus_entropy_d));
max_len=max(cellfun(@length,plus_entropy_d));
plus_entropy_d_mat=nan*ones(size(plus_entropy_d,1),max_len);
for i=1:size(plus_entropy_d,1)
    len=length(plus_entropy_d{i});
    plus_entropy_d_mat(i,1:len)=plus_entropy_d{i};
end

minus_entropy_d=minus_entropy_d(cellfun(@(x) ~isempty(x),minus_entropy_d));
max_len=max(cellfun(@length,minus_entropy_d));
minus_entropy_d_mat=nan*ones(size(minus_entropy_d,1),max_len);
for i=1:size(minus_entropy_d,1)
    len=length(minus_entropy_d{i});
    minus_entropy_d_mat(i,1:len)=minus_entropy_d{i};
end

%%

figure
W=80;
plot(-W:-1, nanmean(minus_entropy_u_mat(:,end-W+1:end)))
hold on
plot(-W:-1,nanmean(plus_entropy_u_mat(:,end-W+1:end)))
% ylim([0 2])

figure
W=80;
for i=1:3
    subplot(3,1,i)
    plot(i:3:W, nanmean(minus_entropy_u_mat(:,i:3:W)))
    hold on
    plot(i:3:W,nanmean(plus_entropy_u_mat(:,i:3:W)))
end

figure
plot(nanmean(minus_entropy_u_mat(:,1:W)))
hold on
plot(nanmean(plus_entropy_u_mat(:,1:W)))
