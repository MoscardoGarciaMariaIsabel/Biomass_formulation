function [enrichment] = GeneEnrichments(GeneList)

% for all possible genes in Recon
load('recon_genes.mat')

GeneList_all = unique(recon_genes.SYMBOL); % 1729 genes in Recon2
M = numel(GeneList_all);

% for GeneList
% GeneList = GeneList_all(randi([1 1175],1,50)); %for testing

N = numel(GeneList);

load('cancer_genes_datamined.mat')
K=[];
x=[];
for i=1:numel(cancer_genes_names)
    source = lower(cancer_genes{i});
    
    K = [K;sum(ismember(lower(GeneList_all), source))];
    x = [x;sum(ismember(lower(GeneList), source))];
end


enrichment = table(cancer_genes_names,num2cell(1-hygecdf(x-1,M,K,N)),num2cell(repmat(M,822,1)),num2cell(K),num2cell(repmat(N,822,1)),num2cell(x),...
    'VariableNames',{'Database','enrichment','Recon_genes','Recon_cancer_genes','Druglist','DrugList_cancer'});

[B, I] = sort(1-hygecdf(x-1,M,K,N),'ascend');

enrichment = enrichment(I,:);
end