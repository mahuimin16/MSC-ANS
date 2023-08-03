 function [Cri,Z_bestloca,Z_best,Q] = criterion(Z,cri)
[S,I]=sort(Z,'descend');
n = size(S,2);
if cri == "in_e"
    for i=1:n
        for k=1:n-1
            Cri(k,i)=-k*(n-k-1)/(n-1)/(n-1)*(1/k*sum(S(1:k,i))-1/(n-k-1)*sum(S((k+1):(n-1),i))).^2;
        end
        [~,Z_bestloca(i)]=min(Cri(:,i));
        Z_best(i)=S(Z_bestloca(i),i);
        Q(1:n,i) =(Z(1:n,i)<Z_best(i));
    end
end
if cri == "diff"
    for i=1:n
        for k=1:n-2
            Cri(k,i)=S(k+1)-S(k,i);
        end
        [~,Z_bestloca(i)]=min(Cri(:,i));
        Z_best(i)=S(Z_bestloca(i),i);
        Q(1:n,i) =(Z(1:n,i)<Z_best(i));
    end
end