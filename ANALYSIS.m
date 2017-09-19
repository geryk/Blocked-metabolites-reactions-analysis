% ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tol=1e-10;

% get model constrained by the growth medium RPMI 1640 
[model_RPMI1640]=getGrowthMediumConstrainedModel(model);

% get model with blocked all reactions producing/consuming any type of 
% folate metabolite
[rxnsFolate,~,L_folateRxn]=findRxnFolate(model);
model_RPMI1640_delAllFolRxns=changeRxnBounds(model_RPMI1640,rxnsFolate,0,'b');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BlOCKED REACTIONS ANALYSIS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% find blocked reactions in native model 
[minFlux,maxFlux] = fastFVA(model_RPMI1640,0,'max','glpk');
BLOCKEDRXNs_RPMI1640=struct;
BLOCKEDRXNs_RPMI1640.minFlux=minFlux;
BLOCKEDRXNs_RPMI1640.maxFlux=maxFlux;
BLOCKEDRXNs_RPMI1640.L_blockedRxns=minFlux>-tol & maxFlux<tol;

% find blocked reactions due to the folate depletion 
[minFlux,maxFlux] = fastFVA(model_RPMI1640_delAllFolRxns,0,'max','glpk');
BLOCKEDRXNs_RPMI1640.minFlux_delAllFolRxns=minFlux;
BLOCKEDRXNs_RPMI1640.maxFlux_delAllFolRxns=maxFlux;
BLOCKEDRXNs_RPMI1640.L_blockedRxns_delAllFolRxns=minFlux>-tol & maxFlux<tol...
    & BLOCKEDRXNs_RPMI1640.L_blockedRxns==0;

% find blocked reactions for every gene deletion in the model_RPMI1640
BLOCKEDRXNs_RPMI1640_allGenes=findBlockedRxns_forEveryGene(model_RPMI1640,...
    allGenes_inRecon22_annot(:,1),BLOCKEDRXNs_RPMI1640);

% compute measures comparing blocked reactions of individual genes with
% blocked reactions associated with folate depletion
BLOCKEDRXNs_RPMI1640_allGenes=addStatistic(BLOCKEDRXNs_RPMI1640_allGenes,...
    BLOCKEDRXNs_RPMI1640,allGenes_inRecon22_annot,allGenes_inRecon22_annot(:,1),model_RPMI1640,'rxn');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BlOCKED METABOLITES ANALYSIS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% find blocked metabolites in native model 
[blockedCpds,L_blockedCpds,maxFlux]=findBlockedMetabolites(model_RPMI1640,model.mets); 
BLOCKEDCPDs_RPMI1640=struct;
BLOCKEDCPDs_RPMI1640.maxFlux=maxFlux;
BLOCKEDCPDs_RPMI1640.blockedCpds=blockedCpds;
BLOCKEDCPDs_RPMI1640.L_blockedCpds=L_blockedCpds;

% find blocked metabolites due to the folate depletion
target_cpds=setdiff(model.mets,BLOCKEDCPDs_RPMI1640.blockedCpds);
[blockedCpds,L_blockedCpds,maxFlux]=findBlockedMetabolites(model_RPMI1640_delAllFolRxns,target_cpds);
BLOCKEDCPDs_RPMI1640.maxFlux_delAllFolRxns=maxFlux;
BLOCKEDCPDs_RPMI1640.blockedCpds_delAllFolRxns=blockedCpds;
BLOCKEDCPDs_RPMI1640.L_blockedCpds_delAllFolRxns=L_blockedCpds;

% find blocked metabolites for every gene deletion in the model_RPMI1640 
BLOCKEDCPDs_RPMI1640_allGenes=findBlockedMetabolites_forEveryGene(model_RPMI1640,...
    allGenes_inRecon22_annot(:,1),BLOCKEDCPDs_RPMI1640.L_blockedCpds);

% compute measures comparing blocked metabolites of individual genes with
% blocked metabolites associated with folate depletion
BLOCKEDCPDs_RPMI1640_allGenes=addStatistic(BLOCKEDCPDs_RPMI1640_allGenes,...
    BLOCKEDCPDs_RPMI1640,allGenes_inRecon22_annot,allGenes_inRecon22_annot(:,1),model_RPMI1640,'cpd');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST DEPENDENCY OF SAM ACCUMULATION ON THE PURINE BIOSYNTHESIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%remove folate dependency of the two biosynthetic steps in the purine
%pathway catalyzed by enzymes: GARFT and AICART
model_RPMI1640_delAllFolRxns_noPurinesFolDep=remove_purinePathwayFolateDependency(model_RPMI1640_delAllFolRxns);

%test if the accumulation of SAM is blocked 
[blockedCpds,~,~]=findBlockedMetabolites(model_RPMI1640_delAllFolRxns_noPurinesFolDep,{'amet[c]'});





