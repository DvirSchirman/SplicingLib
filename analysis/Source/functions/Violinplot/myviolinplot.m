function  h=myviolinplot( data, color_map, d_min, d_max )
%% a wrapper for violinplot that accepts 2D data with different
%% vector size through a cell array

    M=length(data);
    cats=repmat([''],M,1);
    if nargin<2
        if M<6
            color_map=cbrewer('qual','Set1',max([M,3]));
        elseif M<9
            color_map=cbrewer('qual','Set1',M+1);
            color_map(6,:)=[];
        elseif M<11
            color_map=cbrewer('qual','Set1',M+1);
            color_map(7,:)=[];
        else
            color_map=cbrewer('qual','Set1',M+1);
            color_map(8,:)=[];
        end
        
    end
        
    if nargin<4
        d_min=inf;
        d_max=-inf;
        for i=1:M
            if min(data{i})<d_min
                d_min=min(data{i});
            end
            if max(data{i})>d_max
                d_max=max(data{i});
            end
        end
    end
        
%     figure
    hold on
%     h=[];
    data(cellfun(@isempty,data))={randn(2,1)+(d_min-1e4)};
%     data(cellfun(@(x) nansum(x==x(1))==sum(~isnan(x)),data))={[data{cellfun(@(x) nansum(x==x(1))==sum(~isnan(x)),data)}; randn(2,1)+(d_min-1e4)]};
    for i=1:length(data)
        x=data{i};
        if nansum(x==x(1))==sum(~isnan(x))
            data(i)={[data{i}; randn(2,1)+(d_min-1e4)] };
        end
    end
    for i=1:M
        
        v=violinplot([randn(length(data{i}),i-1)+(d_min-1e4),data{i},randn(length(data{i}),M-i)+(d_min-1e4)]);
        v(i).ViolinColor=color_map(i,:);
        v(i).BoxColor=[0 0 0];
%         v(i).MedianColor=[0 0 0]+0.2;
        h(i)=v(i);
    end
    ylim([d_min d_max])
    ax=gca;
    ax.set('XTickLabel',repmat('',M,1))
end

