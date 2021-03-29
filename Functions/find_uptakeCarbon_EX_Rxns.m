function [up_exRxns, up_ex_mets_carbon,uptake_exRxnsInd]=find_uptakeCarbon_EX_Rxns(model, biomass_reaction)
exRxnsInd=[];
uptake_exchange_mets=[];
exRxnsInd=find(sum(abs(model.S),1)==1);
biomass_id= find(ismember(model.rxns, biomass_reaction));
if ~isempty (biomass_id)
   exRxnsInd= setdiff(exRxnsInd,biomass_id);
end

uptake_exRxnsInd = [];
for a = 1:numel(exRxnsInd)
    if ~((model.ub(exRxnsInd(a))&&model.lb(exRxnsInd(a))>0) || (model.ub(exRxnsInd(a))>0&&model.lb(exRxnsInd(a))==0))
        uptake_exRxnsInd = [uptake_exRxnsInd,exRxnsInd(a)];
    end
end
    
for i=1:numel(uptake_exRxnsInd);
    exmets= find(model.S(:,uptake_exRxnsInd(i)));
 
    uptake_exchange_mets(end+1)=exmets;
    if model.S(exmets,uptake_exRxnsInd(i))==1;
        model.S(exmets,uptake_exRxnsInd(i))=-1;
    end
end
up_exRxns=model.rxns(uptake_exRxnsInd);
ex_mets_carbon =  (regexp(model.metFormulas(uptake_exchange_mets),'C'));
is_organic = ~cellfun('isempty', ex_mets_carbon);
up_ex_mets_carbon=model.mets(uptake_exchange_mets(is_organic));


