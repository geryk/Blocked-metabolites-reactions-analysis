function [blockedRxnsStruct]=findBlockedRxns_forEveryGene(model,genes,BLOCKED_other)
tol=1e-10;
nGenes=length(genes);
nRxns=length(model.rxns);
minFluxMatrix=zeros(nRxns,nGenes);
maxFluxMatrix=zeros(nRxns,nGenes);

for i=1:nGenes
    [modelDel,hasEffect,~,~] = deleteModelGenes(model,genes{i}); 
    if hasEffect==true
       [minFlux,maxFlux] = fastFVA(modelDel,0,'max','glpk');
       minFluxMatrix(:,i)=minFlux;
       maxFluxMatrix(:,i)=maxFlux;
    else
       minFluxMatrix(:,i)=BLOCKED_other.minFlux;
       maxFluxMatrix(:,i)=BLOCKED_other.maxFlux;
    end    
end

%find blocked reactions
L_blockedRxns_matrix=minFluxMatrix>-tol & maxFluxMatrix<tol;

%keep information about blocation only for reactions not blocked in native model
L_blockedRxns_matrix(BLOCKED_other.L_blockedRxns,:)=0;

blockedRxnsStruct=struct;
blockedRxnsStruct.targetRxns=model.rxns(~BLOCKED_other.L_blockedRxns);
blockedRxnsStruct.minFluxMatrix=minFluxMatrix;
blockedRxnsStruct.maxFluxMatrix=maxFluxMatrix;
blockedRxnsStruct.L_blockedRxns_matrix=L_blockedRxns_matrix;
end