function [blockedCpdsStruct]=findBlockedMetabolites_forEveryGene(model,genes,L_blockedCpdsDefault)

cpds=setdiff(model.mets,model.mets(L_blockedCpdsDefault));

nCpdsModel=length(model.mets);
nGenes=length(genes);
nCpds=length(cpds);
maxFluxes_matrix=zeros(nCpds,nGenes);
L_blockedCpds_matrix=false(nCpdsModel,nGenes);

parfor i=1:nGenes
    [modelDel,hasEffect,~,~] = deleteModelGenes(model,genes{i});
    
    if hasEffect==true
       [~,L_blockedCpds,fluxes]=findBlockedMetabolites(modelDel,cpds); 
       maxFluxes_matrix(:,i)=fluxes;
       L_blockedCpds_matrix(:,i)=L_blockedCpds;
    end   
end

blockedCpds_matrix=getBlockedCpds(model,L_blockedCpds_matrix,nGenes,nCpds);

blockedCpdsStruct=struct;
blockedCpdsStruct.targetCpds=cpds;
blockedCpdsStruct.maxFluxes_matrix=maxFluxes_matrix;
blockedCpdsStruct.blockedCpds_matrix=blockedCpds_matrix;
blockedCpdsStruct.L_blockedCpds_matrix=L_blockedCpds_matrix;
end

function [blockedCpds_matrix]=getBlockedCpds(model,L_blockedCpds_matrix,nGenes,nCpds)
blockedCpds_matrix=cell(nCpds,nGenes);
for i=1:nGenes
    blockedCpds=model.mets(L_blockedCpds_matrix(:,i));
    blockedCpds_matrix(1:length(blockedCpds),i)=blockedCpds;
end
end
