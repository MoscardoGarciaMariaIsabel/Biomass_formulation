function [enrichment] = GeneEnrichments(GeneList)

% for all possible genes in Recon
load('recon_genes.mat')

GeneList_all = unique(recon_genes.SYMBOL); % 1729 genes in Recon2
M = numel(GeneList_all);

% for GeneList
% GeneList = GeneList_all(randi([1 1175],1,50)); %for testing

N = numel(GeneList);

load('cancer_genes_datamined.mat')
%We want only CRC cell lines
cancer_genes_names_CRC = cancer_genes_names(contains(cancer_genes_names,'colorectal'));
cancer_genes_CRC = cancer_genes(find(contains(cancer_genes_names,'colorectal')));
tissue = {'Colorectal_cell_lines'};

% We want to check if our predicted essential genes (GeneList) are present
% in any of the CRC cell line list. Hence we pull all the CRC cell line EG.

all_genes = [];
for i=1:numel(cancer_genes_names_CRC)
    all_genes = unique([all_genes;lower(cancer_genes_CRC{i})]);
end

K=[];
x=[];    
K = [K;sum(ismember(lower(GeneList_all), all_genes))];
x = [x;sum(ismember(lower(GeneList), all_genes))];


enrichment = table(tissue,num2cell(1-hygecdf(x-1,M,K,N)),num2cell(M),num2cell(K),num2cell(N),num2cell(x),...
    'VariableNames',{'Cell_lines','enrichment','Recon_genes','Recon_cancer_genes','Predicted_EG','Known_EG'});

[B, I] = sort(1-hygecdf(x-1,M,K,N),'ascend');

enrichment = enrichment(I,:);
end