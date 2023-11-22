from setup import *

common_trans_dist = list(set(abc_trans).intersection(abc_h3k27ac_trans, abc_dnase_trans, abc_comb_sqrt_trans))
common_trans_dist # abc_dnase_trans identified two more samples than the other methods

common_trans_corr = list(set(corr_1_trans).intersection(corr_2_trans, corr_3_trans))
common_trans_corr # all of the correlation methods identified the same samples

common_trans_corr = list(set(common_trans_dist).intersection(common_trans_corr))
common_trans_corr # 10 samples identified by both abc distance and correlation methods (results depend on thresholds chosen)