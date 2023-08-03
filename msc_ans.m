function [Z,Z_bestloca,Q,Beta,label,obj] = msc_ans(X,Y,gamma,cri)
% X      : n*di
%% initialize
maxIter = 25 ; % the number of iterations

numclass = length(unique(Y)) ;
numview = length(X);
numsample = size(Y,1);
conD = zeros(numsample,numsample);
D = cell(numview,1);   
Z = zeros(numsample,numsample) ; % n * n
Beta = zeros(numview,1) ;        % v * 1
sigma = 0;
for i = 1:numview
    di = size(X{i},2);
    X{i} = mapminmax(X{i}',0,1); % turn into d*n
    sigma = sigma + 1/numview*sum(sum(cov(X{i})));
    Beta(i) = 1/numview;
    D{i} = pdist(X{i}');
    D{i} = squareform(D{i});
    conD = conD + Beta(i).^2 * D{i};  
end


gamma1 = sigma * gamma;
fconD = exp(-conD/gamma1);
fconD = fconD - diag(diag(fconD));
sumcolfconD = sum(fconD);
sumrowfconD = sum(fconD,2);
for i = 1:numsample
    for j = 1:numsample
        Z(i,j) = 2 * fconD(i,j)/(sumcolfconD(i)+sumrowfconD(j)-2*fconD(i,j));
    end
end

Q = ones(numsample,numsample);
Cri = zeros(numsample,numsample-1);
flag = 1;
iter = 0;

%%
while flag
    iter = iter +1;
    
    
    Z = Z + eye(numsample,numsample);
    term1 = sum(sum(conD .* Z )) + gamma1 * sum(sum(Z.* log(Z))) + sum(min(Cri));
    obj(iter) = term1;
    Z = Z - eye(numsample,numsample);

%     if (iter>1)&&(abs((obj(iter-1)-obj(iter)))<1e-3 || obj(iter)-obj(iter-1)>0|| iter>maxIter || abs(obj(iter)) < 1e-10)
    if (iter>3)&&((obj(iter-1)-obj(iter)>0)&&(obj(iter)==obj(iter-2)) || iter>maxIter  )
        flag = 0;
    end
    
    
    %% optimize beta
    
    for iv = 1:numview
        M(iv) = norm( D{iv}.*Z.*Q , 'fro' );
        M_inv(iv) = 1/M(iv);
    end
    Beta = M_inv ./ (sum(M_inv));
    
    %% optimize q_ij
    [S,I] = sort(Z,'descend');
    [Cri,Z_bestloca,Z_best,Q] = criterion(Z,cri) ; 
    Q = (Q + Q')/2;
    Q = Q < 1;
    
    %% optimize Z
    
    conD = zeros(numsample,numsample);
    for iv = 1:numview
        conD = conD + Beta(iv).^2 .* D{iv}.*Q; 
    end

    fconD = exp(-conD/gamma1);
    fconD = fconD - diag(diag(fconD));
    sumcolfconD = sum(fconD);
    sumrowfconD = sum(fconD,2);
    for iz = 1:numsample
        for jz = 1:numsample
            Z(iz,jz) = 2*fconD(iz,jz)/(sumcolfconD(iz)+sumrowfconD(jz)-2*fconD(iz,jz)) ;
        end
    end
        
    [label] = SpectralClustering(Z,numclass,Y);
    resiter{iter} = Clustering8Measure(Y, label);
           
end
    

 