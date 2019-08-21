#!/bin/bash
# $ -S /bin/bash
# $ -N job-simulation_mCAVIAR_hhuang_1
# $ -cwd = run from this current working directory (relative paths)
# $ -o stdout-Simulation_1_rho.out
# $ -l h_data=6G,h_rt=12:00:00
# $ -t 1-10:1

source ~/.bash_profile
SGE_TASK_ID=1
#rho=(0.95 0.9 0.85 0.8 0.75 0.7 0.65 0.6 0.55 0.5)
#tau=(1 2 3 4 5)
#sigma_g=(5 7 9 11) 
#this_tau=${tau[$SGE_TASK_ID-1]}
#this_sigma_g=${sigma_g[$SGE_TASK_ID-1]}

#simulation='Simulation_result'
#echo $simulation
#mkdir $simulation
#cd $simulation

#outfile='sim_tau_'$SGE_TASK_ID'_results'

#for k in {1..100};do
#for i in {1..100}; do
#python3 caviar.py -l data/ASN_LD/'result_'${SGE_TASK_ID}'simulate_ld.txt' -z data/ASN_Z/'result_'${SGE_TASK_ID}'simulate_z.txt' -o 'result_asn_'${SGE_TASK_ID}
#python3 caviar.py -l data/EUR_LD/'result_'${SGE_TASK_ID}'simulate_ld.txt' -z data/EUR_Z/'result_'${SGE_TASK_ID}'simulate_z.txt' -o 'result_eur_'${SGE_TASK_ID}
#python3 intersect.py -e 'result_eur_'${SGE_TASK_ID}'_set.txt' -a 'result_asn_'${SGE_TASK_ID}'_set.txt' -o 'result_intersect_'${SGE_TASK_ID}'.txt'
g++ Mcaviar.cpp -o Mcaviar
./Mcaviar -l data/EUR_LD/'result_'${SGE_TASK_ID}'simulate_ld.txt' -z data/EUR_Z/'result_'${SGE_TASK_ID}'simulate_z.txt' -a data/ASN_LD/'result_'${SGE_TASK_ID}'simulate_ld.txt' -b data/ASN_Z/'result_'${SGE_TASK_ID}'simulate_z.txt' -o 'result_m_'${SGE_TASK_ID}
		#python3 MCaviar.py -l LD -z Z -o $'results_'$i'_'$j -t $i -s $j
		#python3 MCaviar.py -l data/EUR_LD/'result_'${SGE_TASK_ID}'simulate_ld.txt' -z data/EUR_Z/'result_'${SGE_TASK_ID}'simulate_z.txt' -a data/ASN_LD/'result_'${SGE_TASK_ID}'simulate_ld.txt' -b data/ASN_Z/'result_'${SGE_TASK_ID}'simulate_z.txt' -o 'result_'${i}'_'${j}'_round'${SGE_TASK_ID} -t ${i} -s ${j}
	#python3 simulate_hardcode.py -l ld_out_ra_eur_first50.ld -z ld_out_ra_eur_first50.zscores -o 'sim_z_'$i
#python3 ../MCaviar.py -l LD -z Z -o $outfile'_multiple_'$i -t $this_tau
#python3 ../concatenate.py -z Z
#python3 ../caviar.py -l ../sample_data/50_LD.txt -z concatenated_z_score.txt -o $outfile'_original_'$i
#rm -r LD
#rm -r Z

speculate_file_m='result_m_'${SGE_TASK_ID}'_set.txt'
#speculate_file_o=$outfile'_original_'$i'_set.txt'
python3 capture.py -s $speculate_file_m -t true_causal_set.txt -p 'result_m_'${SGE_TASK_ID}
#python3 ../capture.py -s $speculate_file_o -t true_causal_set.txt -p '../'$simulation'_o'
		#done
#done

#mv '../'$simulation'_m_config_size.txt' '../'Config_Size_m
#mv '../'$simulation'_m_recall_rate.txt' '../'Recall_Rate_m
#mv '../'$simulation'_o_config_size.txt' '../'Config_Size_o
#mv '../'$simulation'_o_recall_rate.txt' '../'Recall_Rate_o