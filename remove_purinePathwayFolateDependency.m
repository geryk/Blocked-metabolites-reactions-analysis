function model=remove_purinePathwayFolateDependency(model)

 %remove GARFT reaction (gene GART) (10fthf[c] + gar[c]  <=> fgam[c] + h[c] + thf[c])
 model=removeRxns(model,'GARFT');
 
 %add new GARFT reaction without dependency on folate (gar[c]  <=> fgam[c] + h[c])
 rxn_def={'gar[c]','<==>','fgam[c]', '+', 'h[c]'};
 rxn_def=strjoin(rxn_def,' ');
 model=addReaction(model,'GARFT',rxn_def);
 
 
 %remove AICART reaction (gene ATIC) (10fthf[c] + aicar[c]  <=> fprica[c] + thf[c] )
 model=removeRxns(model,'AICART');
 
 %add new AICART reaction without dependency on folate (aicar[c] <==> fprica[c])
 rxn_def={'aicar[c]', '<==>', 'fprica[c]'};
 rxn_def=strjoin(rxn_def,' ');
 model=addReaction(model,'AICART',rxn_def);
end