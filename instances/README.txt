You are free to use the instances in this zip file, if you cite the following work:

@article{Levorato2020,
author = {Levorato, M. and Figueiredo, R. and Frota, Y. and Jouglet, A. and Savourey, D.},
journal = {To be published},
keywords = {Microgrid, Energy planning, Robust optimization, Flexible commitments, Real time command strategy},
title = {{Robust energy trading and scheduling for microgrids based on the Contract Collaboration Problem}},
year = {2020}
}


Instances description
=================================
1) full_intances: the 3 microgrid instances used in the main experiments:
  - A_instance2_1NDU_Cons : the *Cons* microgrid instance, with only the consumer uncertain system;
  - A_instance2_1NDU_Prod : the *Prod* microgrid instance, with only the producer uncertain system;
  - A_instance2 : the *Prod\&Cons* microgrid instance, with both uncertain systems (producer and consumer).

  For each instance, there are 3 files, each one for a group of randomly generated scenarios:
  - <instance_file_prefix>_1000s_skewed-left.txt : left-skewed beta distribution with parameters $\alpha=5,\beta=2$;
  - <instance_file_prefix>_1000s_skewed-right.txt : right-skewed beta distribution with parameters $\alpha=2,\beta=5$;
  - <instance_file_prefix>_1000s_uniform.txt: uniform distribution.
  The randomly-generated scenario values are located in the end of each file.

2) toy_sensitivity: reduces microgrid instances used in the sensitivity tests. Each file corresponds to a specific
                     combination of model cost parameters, and the filename follows this mask:
                     OC<OC_cost>_Ct_<fixed_contract_price>_DS<DS_cost>_ST<ST_cost>_NDU_<NDU_profile>_<distribution_type>.txt
                    The parameters OC_cost, DS_cost, ST_cost follow the values listed in Table VI of the aforementioned paper. 
                    (caption: Energy cost parameters in sensitivity analysis). <distribution_type> can be one of the 3 distributions
                    used in the work (skewed-left, skewed-right or uniform). The <fixed_contract_price> used was 'default'.
                    The <NDU_profile> used was 'default'. 
                    Again, the randomly-generated scenario values are located in the end of each file.

