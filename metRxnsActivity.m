function[isBlockedFlux]=metRxnsActivity(BLOCKED_genes,BLOCKED_other,indx_gene,model,cpd)

[rxnList, ~] = findRxnsFromMets(model,cpd);
rxnID = findRxnIDs(model,rxnList);

Lblocked_native=BLOCKED_other.L_blockedRxns(rxnID,:);
Lblocked_gene=BLOCKED_genes.L_blockedRxns_matrix(rxnID,indx_gene);
Lblocked_gene=Lblocked_gene(Lblocked_native==0);

%fractActviveRxns=nnz(Lblocked_gene==0)/length(Lblocked_gene);
if min(Lblocked_native)==1
   isBlockedFlux='blocked in native model';
elseif min(Lblocked_gene)==1
   isBlockedFlux=1;
else   
    isBlockedFlux=0;
end    
end

    

