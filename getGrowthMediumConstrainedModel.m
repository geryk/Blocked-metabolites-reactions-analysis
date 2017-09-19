function [model_medium]=getGrowthMediumConstrainedModel(model)
[~,selUpt] = findExcRxns(model);

%set desired growth medium %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
uptakeRxns=getUptRxns_RPMI1640();
%uptakeRxns=getUptRxns_minimalMedium();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NEW MODEL WITH ALL UPTAKE REACTIONS NOT CONTAINED IN
% GROWTH MEDIUM BLOCKED 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L_medium=ismember(model.rxns,uptakeRxns);
L_Upt_notInMedium=(L_medium==0) & (selUpt==1);
uptakeRxns_notInMedium=model.rxns(L_Upt_notInMedium);
model_medium=changeRxnBounds(model,uptakeRxns_notInMedium,0,'l');
end

function [uptakeRxns]=getUptRxns_RPMI1640()
aminoacids={'EX_gly(e)';'EX_arg_L(e)';'EX_asn_L(e)';'EX_asp_L(e)';'EX_Lcystin(e)';
            'EX_glu_L(e)';'EX_gln_L(e)';'EX_his_L(e)';'EX_4hpro(e)';'EX_ile_L(e)';
            'EX_leu_L(e)';'EX_lys_L(e)';'EX_met_L(e)';'EX_phe_L(e)';'EX_pro_L(e)';
            'EX_ser_L(e)';'EX_thr_L(e)';'EX_trp_L(e)';'EX_tyr_L(e)';'EX_val_L(e)'};
        
vitamins={'EX_btn(e)';'EX_chol(e)';'EX_pnto_R(e)';'EX_ncam(e)';
          'EX_bz(e)';'EX_pydxn(e)';'EX_ribflv(e)';'EX_thm(e)';'EX_adpcbl(e)';
          'EX_inost(e)';'EX_fol(e)'};
      
inorganic={'EX_ca2(e)';'EX_so4(e)';'EX_cl(e)';'EX_na1(e)';'EX_hco3(e)';
            'EX_pi(e)';'EX_k(e)';'EX_o2(e)';'EX_h(e)';'EX_h2o(e)';'EX_nh4(e)';
            'EX_fe3(e)';'EX_fe2(e)';'EX_lnlc(e)';'EX_lnlnca(e)'};
        
other={'EX_glc(e)';'EX_gthrd(e)';'EX_4nph(e)';'EX_lnlc(e)';'EX_lnlnca(e)'};

uptakeRxns=[aminoacids;inorganic;vitamins;other];
end

function [uptakeRxns]=getUptRxns_minimalMedium()
aminoacids={'EX_his_L(e)';'EX_ile_L(e)';'EX_lys_L(e)';'EX_met_L(e)';'EX_thr_L(e)';'EX_trp_L(e)';'EX_val_L(e)'};
essentialKinetensinAmk={'EX_pro_L(e)';'EX_tyr_L(e)';'EX_leu_L(e)'};
inorganic={'EX_pi(e)'};

uptakeRxns=[aminoacids; essentialKinetensinAmk; inorganic];
end