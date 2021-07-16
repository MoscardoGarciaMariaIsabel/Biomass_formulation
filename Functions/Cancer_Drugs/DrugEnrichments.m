function [enrichment,A] = DrugEnrichments(DrugList)

% for all possible drugs in Recon
% load('P:\PhD Project\Resources\DrugBank\v5.1.3\dico_drugs_final_recon_filtered.mat')
load dico_drugs_final_recon_filtered.mat


dico_drugs = dico_final_recon_filtered;

DrugList_all = unique(dico_drugs.DrugName(contains(dico_drugs.Action,'inhib'))); % 1175 drugs to test
DrugList_all = lower(DrugList_all);
M = numel(DrugList_all);
% DrugList_all = regexprep(DrugList_all,'\W','')
% for DrugList
% DrugList = DrugList_all(randi([1 1175],1,50)) %for testing

DrugList = lower(DrugList);
N = numel(DrugList);
% DrugList = regexprep(DrugList,'\W','')

load('cancer_drugs_datamined.mat')
K=[];
x=[];

A = zeros(numel(DrugList),numel(cancer_drugs_names));
for i=1:numel(cancer_drugs_names)
    source = lower(cancer_drugs{i});
    
    %split if there are parentheses
    source2=[];
    for j=1:numel(source)
        source2 = [source2;strsplit(source{j},'(')'];
    end
    source2 = regexprep(source2,'\)|Hydrochloride|Liposome','');
    source2 = strtrim(source2);
    source2 = unique(source2);
%     source2 = regexprep(source2,'\W','')
    K = [K;sum(ismember(DrugList_all, source2))];
    x = [x;sum(ismember(DrugList, source2))];
    
    A(ismember(DrugList, source2),i) = 1;
end
sum(A,2);

enrichment = table(cancer_drugs_names,num2cell(1-hygecdf(x,M,K,N)),num2cell(repmat(M,14,1)),num2cell(K),num2cell(repmat(N,14,1)),num2cell(x),...
    'VariableNames',{'Database_website','enrichment','Recon_drugs','Recon_cancer_drugs','Druglist','DrugList_cancer'});



end