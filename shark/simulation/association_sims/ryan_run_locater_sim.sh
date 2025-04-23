#! /bin/bash

# Resources Required
# for locater sims:
# 80Gb and 4 threads needed with no checkpoints (as in locater 11 and 12)
# 160Gb and 8 threads needed with 3 checkpoints (as in locater 8, locater 10, locater 13, locater14, locater16)
# 10Gb and 1 thread required for staar runs (staar14-staar17) -- got some through with 5Gb
# 5 Gb and 1 thread required for ARG-NEEDLE
mem=160
nthreads=8

declare -a run_array=("locater17")

# image="docker(rchrist7/mini-shark-no-staar)" # used for locater runs
# image="docker(rchrist7/gdistill)" # used for staar runs
image="docker(rchrist7/mini-shark)" # used for argneedle runs

for run_name in "${run_array[@]}";
do

run_group="locater_sims"

r_script="/rchrist/shark/simulation/association_sims/ryan_sim_parameters/${run_name}.R"

group="/rchrist/${run_group}/${run_name}"

bgadd "$group"

command="/storage1/fs1/ccdg/Active/analysis/rchrist_test/shark/simulation/association_sims/locater_sims"

queue="ccdg"

logdir="/storage1/fs1/ccdg/Active/analysis/rchrist_test/logs/${run_group}/${run_name}"

#already_done=(1 232 233 234 235 236 237 238 239 240 241 242 243 244 245 246 247 248 419 870 872 873 874 875 876 877 878 879 880 881 882 883 884 885 886 887 888 890 891 892 893 894 895 896 897 899 900)


######### BELOW HERE SHOULDN'T CHANGE OFTEN #########

rm -r "$logdir"; mkdir -p "$logdir"


user_group="compute-ccdg"

for iter in {381,487,579,595};
#{1..900};
do

# if [[ $(echo ${already_done[@]} | fgrep -w $iter) ]]
# then
#     continue
# fi

if test -f /storage1/fs1/ccdg/Active/analysis/rchrist_test/data/simulation/locater_sims/$run_name/res_$iter.rds; then
    continue
fi


jobname="${run_name}_$iter"

selectString="select[mem>${mem}GB] rusage[mem=${mem}GB] span[hosts=1]" # affinity[core(${nthreads}):membind=localprefer]"
#affinity[core(${nthreads},same=numa):cpubind=numa:membind=localonly]" membind=localprefer

export HOME=/home/${USER}
export LSF_DOCKER_VOLUMES="/storage1/fs1/ccdg/Active/analysis/rchrist_test:/rchrist /home/rchrist:/home/rchrist"
export DISPLAY=${DISPLAY}
# -G flag used to be ACCOUNTS//compute-ccdg
#bsub -oo "$logdir/log_$iter.out" -m "$host_string" -g "$group" -G "$user_group" -q "$queue" -J "$jobname" -M "${mem}GB" -R "$selectString" -a "$image" /bin/bash -x "$command" "$nthreads" "$iter"
bsub -oo "$logdir/log_$iter.out" -n "$nthreads" -g "$group" -G "$user_group" -q "$queue" -J "$jobname" -M "${mem}GB" -R "$selectString" -a "$image" /bin/bash -x "$command" "$r_script" "$nthreads" "$iter"

echo $iter

sleep 0.1

done
done

