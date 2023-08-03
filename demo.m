clear;
clc;
warning off;
addpath(genpath('./'));

%% dataset
% ds={'Yale','yaleA_3view','3sources','Texas','Cornell','WebKB_cor2views','MSRCV1','Washington','WebKB_Wisconsin2views'};
% ds={'Wisconsin','Dermatology','ORLRnSp','ORL_3Views','ORL_4Views','NGs','ForestTypes'};
% ds={'BBCSport','prokaryotic','synthetic3d','Movies','BBC','WikipediaArticles','proteinFold'};
% ds={'WebKB','Reuters-1200','buaaRnSp','Flower17','Caltech101-7','100Leaves'};
% ds={'HW_2Views','HW_6Views','MFeat_2Views','uci-digit','Caltech101-20','BDGP','Cora','Wiki_fea'};
% ds={'CiteSeer','NUS-WIDE-Lite','NUS-WIDE-SCENE'};
% ds={'BDGP_2V','Caltech-5V','Caltech20_6V','LandUse-21_3V','RGB-D_2V','Scene-15_3V'};
% ds={'MNIST_USPS_2V','NUS-WIDE-OBJECT-10','Reuters-7200','Caltech101','Hdigit','SUNRGBD'};
% ds={'NUS-WIDE-OBJECT-10','Reuters-7200','Caltech101','Hdigit','SUNRGBD'};
% ds={'Hdigit'};
% ds={'MSRCV1','Dermatology','ORLRnSp','ForestTypes','WikipediaArticles','MFeat_2Views','uci-digit','100Leaves','HW_6Views','Wiki_fea','Reuters-7200'};
% ds={'MSRCV1','Dermatology','MFeat_2Views','Caltech101-20','BDGP','Wiki_fea','Reuters-7200'};
ds={'Dermatology'};

% ds={'MSRCV1','Dermatology','ORLRnSp','ForestTypes','WikipediaArticles','MFeat_2Views','uci-digit','100Leaves','HW_6Views','Wiki_fea','Reuters-7200'};


% ds = {'Yale'};
% dsPath = '.\datasets\';
dsPath = 'D:\Master\Research\Code\Multi-view Data\';
resultdir = '.\res\';
metric = {'ACC','nmi','Purity','Fscore','Precision','Recall','AR','Entropy'};
% gamma = 1;
gamma = [1/20];
cri = ["in_e"];

% cri = ["in_e","diff"];


for dsi =1:length(ds)
    dataName = ds{dsi}; disp(dataName);
    load(strcat(dsPath,dataName));
%     X= data';
%     Y= truth ;
    k = length( unique(Y));
    n = length(Y);
    %%
    for dcri = 1:length(cri)
    for id = 1:length(gamma)       
        tic;
        [Z,Z_bestloca,Q,Beta,label,obj] = msc_ans(X,Y,gamma(id),cri(dcri));
        res = Clustering8Measure(Y, label);
        betaall{id,dcri} = Beta;
        locaall{id,dcri} = Z_bestloca;
        timer(id,dcri)  = toc;
        resall{id,dcri} = res;
        objall{id,dcri} = obj;
    end
    end
%     save([resultdir, char(dataName),'_result.mat'], 'resall', 'objall','betaall','timer','locaall');
end


