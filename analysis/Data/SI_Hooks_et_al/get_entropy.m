function e = get_entropy(msa)
    
    Alphabet = 'ACGT';
    
    if ~iscell(msa)
        msa=mat2cell(msa,ones(size(msa,1),1));
    end
    msa = cellfun(@upper,msa,'un',0);
    min_len=min(cellfun(@length,msa));
    msa = cellfun(@(x) x(1:min_len),msa,'un',0);
    
    prob_mat=zeros(4,min_len);
    
    for i=1:size(prob_mat,2)
        for j=1:4
            prob_mat(j,i)=sum(cellfun(@(x) x(i)==Alphabet(j),msa))/size(msa,1);
            a=1;
        end
    end
    e=-nansum(prob_mat.*log2(prob_mat));
    


end

