function BLOCKED_genes=addStatistic(BLOCKED_genes,BLOCKED_other,allGenesInModel_annot,genes,model,type)

%Set variables according to type of analysis
if strcmp(type,'cpd')
   L_blockedMatrix=BLOCKED_genes.L_blockedCpds_matrix;
   L_blocked_delAllFolRxns=BLOCKED_other.L_blockedCpds_delAllFolRxns;
elseif strcmp(type,'rxn')
   L_blockedMatrix=BLOCKED_genes.L_blockedRxns_matrix;
   L_blocked_delAllFolRxns=BLOCKED_other.L_blockedRxns_delAllFolRxns;
end    
   
Lnonzero=sum(L_blockedMatrix,1)>0;
IndxNonzero=find(Lnonzero);
nNonzero=length(IndxNonzero);
nBlocked_delAllFolRxns=nnz(L_blocked_delAllFolRxns);

Stat=table;
Stat.HGNC_id=genes(Lnonzero);
gene=cell(nNonzero,1);
enzName=cell(nNonzero,1);
nBlocked_nonzero=zeros(nNonzero,1);

%Measures comparing blocked metabolites associted with gene deltions with 
%set of blocked metabolites associated with folate depletion
IntersectFract_gene=zeros(nNonzero,1);
IntersectFract_folDepl=zeros(nNonzero,1);
Jaccard=zeros(nNonzero,1);

%Measures comparing blocked metabolites associted with gene deltions with 
%metabolites essential for main folate functions: neurotransmiter synthesis,
%methylation reactions.
targets={'thbpt[c]';'amet[c]'};
colNames={'BH4_blockAccumulation';'BH4_blockTurnover';'SAM_blockAccumulation';'SAM_blockTurnover'};
isBlockedTable=cell(nNonzero,2*length(targets));

for i=1:nNonzero
    i
    %annotation of gene ID with gene and enzyme name 
    genes_i=Stat.HGNC_id{i};
    if isempty(strfind(genes_i,'HGNC'))==1
       [gene_i,enzName_i]=annotGeneNCBI(genes_i,allGenesInModel_annot);
    elseif isempty(strfind(genes_i,'HGNC'))==0
           [gene_i,enzName_i]=annotGeneHGNC(genes_i,allGenesInModel_annot);
    end       
    gene{i}=gene_i;
    enzName{i}=enzName_i;
    
    Lblocked_i=L_blockedMatrix(:,IndxNonzero(i));
    nBlocked_nonzero(i,1)=nnz(Lblocked_i);
    
    %Measures comparing blocked metabolites associted with gene deletions with 
    %set of blocked metabolites associated with folate depletion
    nIntersect=nnz(Lblocked_i & L_blocked_delAllFolRxns);
    nUnion=nnz(Lblocked_i | L_blocked_delAllFolRxns);
    IntersectFract_gene(i,1)=nIntersect./nBlocked_nonzero(i);
    IntersectFract_folDepl(i,1)=nIntersect./nBlocked_delAllFolRxns;
    Jaccard(i,1)=nIntersect./nUnion;
    
    %Test if the metabolites essential for main folate functions are
    %blocked in the case of genes_i deletion
    isBlockedTable=getIsBlockedTable(targets,isBlockedTable,i,IndxNonzero(i),model,BLOCKED_genes,BLOCKED_other,type);   
end

Stat.gene=gene;
Stat.enzName=enzName;
Stat.nBlocked=nBlocked_nonzero;
Stat.IntersectFract_gene=IntersectFract_gene;
Stat.IntersectFract_folDepl=IntersectFract_folDepl;
Stat.Jaccard=Jaccard;
isBlockedTable=cell2table(isBlockedTable, 'VariableNames',colNames);
Stat=[Stat isBlockedTable];   

BLOCKED_genes.Stat=Stat;
end

function [isBlocked,isBlockedFlux]=blockInfo_cpd(model,BLOCKED_genes,BLOCKED_other,indx_col,cpd,type)

%test if the flux through the metabolite "cpd" is blocked
if strcmp(type,'rxn')
   isBlockedFlux=metRxnsActivity(BLOCKED_genes,BLOCKED_other,indx_col,model,cpd);
   isBlocked='-';
   
   %test if production of the metabolite "cpd" is blocked
elseif strcmp(type,'cpd');
   indx_BH4=find(ismember(BLOCKED_genes.targetCpds,cpd));
   indx_BH4_model=findMetIDs(model,cpd); 
   isBlockedFlux='-';
    if isempty(indx_BH4)==0
       isBlocked=BLOCKED_genes.L_blockedCpds_matrix(indx_BH4_model,indx_col);
    else
        isBlocked='blocked in native model';
    end
end
end

function [isBlockedTable]=getIsBlockedTable(metsTargets,isBlockedTable,indx_row,indx_col,model,BLOCKED_genes,BLOCKED_other,type)
nMets=length(metsTargets);
for i=1:nMets
    [isBlocked,isBlockedFlux]=blockInfo_cpd(model,BLOCKED_genes,BLOCKED_other,indx_col,metsTargets{i},type);
    isBlockedTable{indx_row,i*2-1}=isBlocked;
    isBlockedTable{indx_row,i*2}=isBlockedFlux;
end
end

%helper function used by two following annotation functions
function [gene_i,enzName_i]=annotFromTable(genes_i,allGenesInModel_annot,colIndx_enzname,colIndx_geneName,colIndx_symb)
genes_i=unique(genes_i);
gene_i='';
enzName_i='';
for j=1:length(genes_i)
    if j>1
        sep=';';
    else
        sep='';
    end    
    L=ismember(allGenesInModel_annot(:,colIndx_symb),genes_i{j});
    gene_i=strcat(gene_i,sep,allGenesInModel_annot{L,colIndx_geneName});
    enzName_i=strcat(enzName_i,sep,allGenesInModel_annot{L,colIndx_enzname});
end
end

%annotate HGNC genes by enz. names
function [gene_i,enzName_i]=annotGeneHGNC(genes_i,allGenesInModel_annot)
genes_i=strsplit(genes_i,';');
[gene_i,enzName_i]=annotFromTable(genes_i,allGenesInModel_annot,4,3,1);
end

%annotate NCBI gene codes by HGNC symbols and enz. names 
function [gene_i,enzName_i]=annotGeneNCBI(genes_i,allGenesInModel_annot)
genes_i=strsplit(genes_i,';');
for j=1:length(genes_i)
    x=strsplit(genes_i{j},'.');
    genes_i{j}=x{1,1};
end
[gene_i,enzName_i]=annotFromTable(genes_i,allGenesInModel_annot,3,2,4);
end


    




