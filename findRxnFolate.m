%find all reactions having any folate metabolite as a substrate or product
function [rxnsFolate,rxnsFolateNames,L_folateRxn]=findRxnFolate(model)
occurence=strfind(lower(model.metNames),'fol');
L_folateCpds=~cellfun(@isempty,occurence);

folateCpds=model.mets(L_folateCpds);

[rxnsFolate, ~] = findRxnsFromMets(model,folateCpds);
L_folateRxn=ismember(model.rxns,rxnsFolate);
rxnsFolateNames=model.rxnNames(L_folateRxn);
end
