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
1) utc_skew: the 3 UTC microgrid instances used in the main experiments:
  - A_instance2_1NDU_Cons : the *Cons* microgrid instance, with only the consumer uncertain system;
  - A_instance2_1NDU_Prod : the *Prod* microgrid instance, with only the producer uncertain system;
  - A_instance2 : the *Prod\&Cons* microgrid instance, with both uncertain systems (producer and consumer).

  For each instance, there are 3 files, each one for a group of randomly generated scenarios:
  - <instance_file_prefix>_1000s_skewed-left.txt : left-skewed beta distribution with parameters $\alpha=5,\beta=2$;
  - <instance_file_prefix>_1000s_skewed-right.txt : right-skewed beta distribution with parameters $\alpha=2,\beta=5$;
  - <instance_file_prefix>_1000s_uniform.txt: uniform distribution.
  The randomly-generated scenario values are located in the end of each file.

2) toy: reduces microgrid instances used in the sensitivity tests. Each file corresponds to a specific
                     combination of model cost parameters, and the filename follows this mask:
                     OC<OC_cost>_Ct_<fixed_contract_price>_DS<DS_cost>_ST<ST_cost>_NDU_<NDU_profile>_<distribution_type>.txt
                    The parameters OC_cost, DS_cost, ST_cost follow the values listed in Table VI of the aforementioned paper. 
                    (caption: Energy cost parameters in sensitivity analysis). <distribution_type> can be one of the 3 distributions
                    used in the work (skewed-left, skewed-right or uniform). The <fixed_contract_price> used was 'default'.
                    The <NDU_profile> used was 'default'. 
                    Again, the randomly-generated scenario values are located in the end of each file.

3) japan_microgrid: Multiyear microgrid data from a research building in Tsukuba, Japan. RCCP problem instances and scenarios were assembled based on the dataset provided at the following website: https://springernature.figshare.com/articles/dataset/Multiyear_microgrid_data_from_a_research_building_in_Tsukuba_Japan/7403954?backTo=/collections/Multiyear_microgrid_data_from_a_research_building_in_Tsukuba_Japan/4176698
*** Reference:
Vink, Karina; Ankyu, Eriko; Koyama, Michihisa (2019): Multiyear microgrid data from a research building in Tsukuba, Japan. figshare. Collection. https://doi.org/10.6084/m9.figshare.c.4176698.v1 
