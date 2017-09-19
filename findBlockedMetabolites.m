function[blockedCpds,L_blockedCpds,fluxes]=findBlockedMetabolites(model,cpds)
%model = changeRxnBounds(model,rxnsDel,0,'b');
changeCobraSolver('glpk');
tol=1e-10;

Lsink=(sum(model.S(:,:),1)==-1 & sum(model.S(:,:)~=0,1)==1);

nCpdsModel=length(model.mets);
nCpds=length(cpds);

L_blockedCpds=false(nCpdsModel,1);
fluxes=zeros(nCpds,1);

for i=1:nCpds
    %i
    cpd_i=cpds{i};
    indx_cpd_i= findMetIDs(model,cpd_i);
    indxRxnSink=find(Lsink & model.S(indx_cpd_i,:)==-1);
    
    if length(indxRxnSink)>1
        indxRxnSink
       sprintf('>1 identical sink reactions found !!!')
       break 
    end
    
    if isempty(indxRxnSink)
       sinkReaction_def=sprintf('%s --> ',cpd_i);
       target_name='targetRxn';
       model=addReaction(model,target_name,sinkReaction_def);
    else
        target_name=model.rxns(indxRxnSink);
    end
    
    model=changeObjective(model,target_name);
    FBAsolution_targetRxn=optimizeCbModel(model, 'max');
    fluxes(i)=FBAsolution_targetRxn.f;
    
    if FBAsolution_targetRxn.f<tol
       L_blockedCpds(indx_cpd_i)=true;
    end
    
    if isempty(indxRxnSink)
       model=removeRxns(model,target_name);
    end   
end  

blockedCpds=model.mets(L_blockedCpds);
end

   
    
    


