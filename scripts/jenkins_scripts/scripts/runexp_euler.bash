#! /bin/bash

# ICON
#
# ---------------------------------------------------------------
# Copyright (C) 2004-2026, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# Contact information: icon-model.org
#
# See AUTHORS.TXT for a list of authors
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ---------------------------------------------------------------

#------------------------------------------------------------------------------
# This script is called on all buildbot slaves and should be used to run the runscripts
#------------------------------------------------------------------------------

# Define of the text which is written in case of a script fail or OK run

TEXT_RUN_FAILED="FAILED"
TEXT_RUN_OK="OK"

#==============================================================================
stop_on_error()
{
# Check if the first parameter (return status) is not OK
    echo STATUS_IN_FILE ${STATUS_IN_FILE}
    if [[ $1 -ne 0 || ${STATUS_IN_FILE} -ne 0 ]]
    then
      if [[ $1 -ne 0 ]]
      then
        printf '%-50s : %s \n' $2 ${TEXT_RUN_FAILED} >> ../${LOOP_STATUS_FILE}
        exit $1
      else
        printf '%-50s : %s\n' $2 ${TEXT_RUN_FAILED} >> ../${LOOP_STATUS_FILE}
        exit ${STATUS_IN_FILE}
      fi
    fi

    printf '%-50s : %s\n' $2 ${TEXT_RUN_OK} >> ../${LOOP_STATUS_FILE}
}

#------------------------------------------------------------------------------------

warning_on_error()
{
# Check if the first parameter (return status) is not OK
    if [[ $1 -ne 0 || ${STATUS_IN_FILE} -ne 0 ]]
    then
      if [[ $1 -ne 0 ]]
      then
        printf '%-50s : %s \n' $2 ${TEXT_RUN_FAILED} >> ../${LOOP_STATUS_FILE}
        echo "*********** WARNING: script failed  $1 ****************"
      else
        printf '%-50s : %s \n' $2 ${TEXT_RUN_FAILED} >> ../${LOOP_STATUS_FILE}
        echo "*********** WARNING: script failed  ${STATUS_IN_FILE} ****************"
      fi
    else
      printf '%-50s : %s\n' $2 ${TEXT_RUN_OK} >> ../${LOOP_STATUS_FILE}
      if [ "${build_post_file}" == "true" ]
      then
        name="`echo $2 | cut -d '.' -f2`"
        RUN_POST="${RUN_POST} post.${name}.run"
        echo "${RUN_POST}" > ./run_post_list

        compare_run="post.${name}_compare_restarts.run"
        if [  -a ./${compare_run} ] ; then
          RUN_POST_COMP="${RUN_POST_COMP} $compare_run"
          echo "${RUN_POST_COMP}" > ./run_post_comp_list
        fi
      fi
    fi
}

#------------------------------------------------------------------------------------

return_ok()
{
    # Exit with status = 0 = OK
    # Arguments:
    #   $1 =  message
    echo "return_ok()"
    echo "$1"

    exit 0
}

#------------------------------------------------------------------------------------



#==============================================================================
# runs the scrpits
run_scripts_submit()
{

  echo "|======================================================|"
  echo "|                                                      |"
  echo "|         Running exp-scripts                          |"
  echo "|                                                      |"
  echo "|======================================================|"

  LOOP_STATUS_FILE="LOOP_STATUS_EXP_FILE"

  rm -f ${LOOP_STATUS_FILE}

  stop_check="warning_on_error"

  echo "Run all *.run in run directory"
  cd run
  EXP_FILES=`cat runscript_list`

  case $target in
      euler*)
           batch_system=slurm
           slurm_user=`whoami`
           typeset -A job_submitted
           if [[ $JENKINS_NO_OF_CORES -lt 13 ]]; then
               memory_per_core=$(expr 12288 / $JENKINS_NO_OF_CORES)
           else
               memory_per_core=1024
           fi
           case $JENKINS_NODE_TYPE in
               6) euler_node=EPYC_7742;;
               7) euler_node=EPYC_7H12;;
           esac

           submmit=""
           submit_details=""
      ;;
      *)
          batch_system=other
      ;;
  esac

  for EXP_FILE in $EXP_FILES
  do
    run_command="$submit $submit_details ./$EXP_FILE"

    if [ -r $EXP_FILE ]
    then
        echo "---------------------------------------------------------"
        echo " Submit new Script: ${EXP_FILE} at $(date)"
        echo " $run_command "
        echo "---------------------------------------------------------"

        [[ ! -f submit.$EXP_FILE ]] ||  rm submit.$EXP_FILE
        echo $run_command >> submit.$EXP_FILE
        if [ "_$submit" = "_qsub -Wblock=true" ]
        then
            echo "echo \$? > ${EXP_FILE}.status.2" >> submit.$EXP_FILE
        fi
        chmod +x submit.$EXP_FILE
        if [[ "$batch_system" = "bsub_euler" ]]
        then
          jobid=$(./submit.$EXP_FILE 2>&1 | sed '/^$/d' | grep "<" | awk '{print $2}'| tr -d "<>")
        else
          ./submit.$EXP_FILE &

    fi


    fi
    if [[ "$batch_system" = "slurm" ]]
    then
        sleep 2
        jobid=`squeue -u ${slurm_user} -h -o '%i' -S '-i' | awk 'NR==1{print $1}'`
        job_submitted["$EXP_FILE"]=$jobid
    fi

    if [[ "$batch_system" = "bsub_euler" ]]
    then
        sleep 2
        echo Watch job with ID $jobid
        job_submitted["$EXP_FILE"]=$jobid
        job_state="RUN"
        while [[ "${job_state}" = "RUN" || "${job_state}" = "WAIT" || "${job_state}" = "PEND" ]]
          do
            job_state=$(bjobs $jobid  | awk 'NR==2{print $3}')
          done

    fi
  done

  # wait for all jobs to finish
  wait
  sleep 30

  echo ${pwd}


  # print and check the results
  for EXP_FILE in $EXP_FILES
  do

    if [ -r $EXP_FILE ]
    then

      # JJ: this part enables the use of out-of-build source for testing, but fails
      # especially on Jenkins

      # Possible fix: Adapt path of final.status file written at the end of file
      # The file is currently written to root/run and not icon_dir/run

      #if [[ "$batch_system" = "slurm" ]]
      #then
      #    sleep 15
      #    slurm_jobid=${job_submitted["$EXP_FILE"]}
      #    job_state=COMPLETING
      #    while [[ "${job_state}" = "COMPLETING" || "${job_state}" = "RUNNING" ]]
      #    do
      #        sleep 5
      #        job_state=`sacct -j ${slurm_jobid} -n -o jobid,state,exitcode | awk 'NR==1{print $2}'`
      #    done
      #    case "$job_state" in
      #        TIMEOUT)
      #            echo 127 > ${EXP_FILE}.final_status
      #        ;;
      #        COMPLETED)
      #            echo 0 > ${EXP_FILE}.final_status
      #        ;;
      #        *)
      #            echo 191 > ${EXP_FILE}.final_status
      #        ;;
      #    esac
      #fi

      exp_status_file=$EXP_FILE.final_status

      if [[ "$batch_system" = "moab_kit" ]]
      then
          sleep 15
          job_state="Running"
          while [[ "${job_state}" = "Idle" || "${job_state}" = "Running" ]]
          do
            job_state=`checkjob $jobid  | awk 'NR==4{print $2}'`
          done
          case "$job_state" in
          Completed*)
                  echo 0 > ${exp_status_file}
              ;;
          esac

      fi

      STATUS_IN_FILE=255
      if [ -r ${exp_status_file} ]
      then
        STATUS_IN_FILE=`cat ${exp_status_file}`
      else
        if [ -r $EXP_FILE.status.2 ]
        then
          STATUS_IN_FILE=`cat ${EXP_FILE}.status.2`
        fi
      fi

      $stop_check 0 $EXP_FILE
    else
      echo "---------------------------------------------------------"
      echo " ${EXP_FILE} does not exist"
      echo "---------------------------------------------------------"
      echo
    fi
  done

  cd ..

  ALL_RUNS_OK=`grep ${TEXT_RUN_FAILED} ${LOOP_STATUS_FILE}`

  if [ $? == 0 ]
  then
    echo "One or more Exp-Runs were not successful"
    RETURN_STATUS=1
  else
    echo "All Exp-Runs were successful"
    RETURN_STATUS=0
  fi



  #==============================================================
  echo "|======================================================|"
  echo "|                                                      |"
  echo "|           Ends                                       |"
  echo "|              $(date)                                 |"
  echo "|                                                      |"
  echo "|======================================================|"
}
#==============================================================================

#------------------------------------------------------------------------------\
# read set-up info
. ./run/set-up.info
submit=$use_sync_submit
target=$use_target
#-----------------------------------------------------------------------------
# load ../setting if exists
if [ -a ./setting ]
then
  . ./setting
fi
#-----------------------------------------------------------------------------
run_scripts_submit
#------------------------------------------------------------------------------
# return OK Status
exit ${RETURN_STATUS}
