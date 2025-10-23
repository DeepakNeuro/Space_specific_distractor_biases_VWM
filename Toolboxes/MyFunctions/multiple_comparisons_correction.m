function [p,sigpvals,nonsigpvals,sigpvalsind,nosigpvalsind] = multiple_comparisons_correction(p,type)

alpha = 0.05;
m = length(p);

%% Holm
if strcmp(type,'holm');
    [psort,pind] = sort(p);
    for i1 = 1:m;
        if psort(i1) > alpha/(m+1-i1);
            R = i1;
            break;
        else
            continue;
        end
    end
    if R == 1,
        fprintf('No hypotheses rejected\n');
    else
        sigpvalsind = pind(1:R-1);
        sigpvals = p(sigpvalsind);
        nosigpvalsind = pind(R:end);
        nonsigpvals = p(nosigpvalsind);
    end
    
end

%% Hochberg
if strcmp(type,'hochberg');
    [psort,pind] = sort(p);
    for i1 = m:-1:1;
        if psort(i1) <= alpha/(m+1-i1);
            R = i1;
            break;
        else
            continue;
        end
    end
    sigpvalsind = pind(1:R-1);
    sigpvals = p(sigpvalsind);
    nosigpvalsind = pind(R:end);
    nonsigpvals = p(nosigpvalsind);
end

%% Benjamni - Hochberg's
if strcmp(type,'Benjaminihochberg');
    [psort,pind] = sort(p);
    for i1 = m:-1:1;
        if psort(i1) <= (i1*alpha)/m;
            R = i1;
            break;
        else
            continue;
        end
    end
    if exist('R');
        sigpvalsind = pind(1:R);
        sigpvals = p(sigpvalsind);
        nosigpvalsind = pind(R+1:end);
        nonsigpvals = p(nosigpvalsind);
    else
        fprintf('No significance! \n');
        sigpvalsind = [];
        sigpvals = [];
        nosigpvalsind = [];
        nonsigpvals = [];
        return;
    end
end


%% Plot

pbin = zeros(m,1);
pbin(sigpvalsind)=1;

figure;
for i1=1:m;
    if pbin(i1)==1;
        line([i1,i1],[0,1],'Color','r');
        s = sprintf('%1.4f',p(i1));
        h1 = text(i1,1.05,s);
        set(h1,'Rotation',90,'FontSize',6);
        hold on;
    end
end
%set (gca,'YTick',1:6,'YTickLabel',{'0','0.2','0.4','0.6','0.8','1'});
ylim([-0.2,1.45]);
line([1,m],[0,0],'Color','r');
s = sort(sigpvalsind);
for j1=1:length(s);
    h2 = text(s(j1),1.26,sprintf('(%d)',s(j1)));
    set(h2,'Rotation',90,'FontSize',6,'FontWeight','bold');
end
xlim([0 m+1]);
xlabel('Hypotheses');
ylabel('Significance');
if strcmp(type,'holm');
    title('Holm Correction')
elseif strcmp(type,'hochberg');
    title('Hochberg Correction');
elseif strcmp(type,'Benjaminihochberg');
    title('Benjamini - Hochberg Correction');
else
end

box on;
end