#!/usr/bin/env python

import os
#import numpy as np
import datetime 
import time
import logging
import math

config = 7 

if config==7:
# paros: motion_ut1
  #dir_parent = "/home/plyu/hosmnt/projects/HornsRev"
  dir_pre = "motion_ut2_"
  pbs_pre = "mo_u2_"
  #fn_log = '/home/plyu/batch_pbs'+'_motion_ut1'+'.log'
  #pbs_script = "run_second_paros.sh"
  #pbs_email = " -M plyu@umn.edu "
  
  dir_parent = "/work/zengx372/Pin/HornsRev/motion_ut2"
  fn_log = '/u/zengx372/Pin/batch_pbs'+'_motion_ut2'+'.log'
  pbs_script = "run_second_copper.sh"
  icase_s = 141
  icase_e = 160
  ncase = icase_e - icase_s + 1
  pbs_email = " "

  ## Motion_ut1, Aegean, 129 z grid case, 2500 steps, 20.3H
  ## Motion_ut1, Paros, 65 grid case, 10m 150 steps, 5h20m 5000 steps
  sleep_first = 60 * 100
  sleep_later = 60 * 10

logging.basicConfig(filename=fn_log, level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')
logging.info('---------------------------------------------------------------------')
logging.info('Batch_PBS: a script to submit tasks automatically after previous task is done')
logging.info('---------------------------------------------------------------------')
logging.info('Specifically, this script run job '+dir_pre+str(icase_s)+" to "+dir_pre+str(icase_e))

t_begin = datetime.datetime.now()
t_tmp = ncase * sleep_first
t_end_est = t_begin + datetime.timedelta(seconds=t_tmp)
logging.info('Finish-time estimated as: '+str(t_end_est)+'\n')


def check_qstat(id_in):
  #result_qstat = os.system("qstat "+str(id_in))
  #print('result_qstat = '+ str(result_qstat)) 
  #if result_qstat == 0:
  #  qstat_out = 1
  #  print("check_qstat: Task "+str(id_in)+" is alive.")
  #else:
  #  qstat_out = 0
  #  print("check_qstat: Task "+str(id_in)+" is dead.")
  
  result_qstat = os.popen("qstat "+str(id_in)).readlines()
  logging.debug('result_qstat =\n '+ str(result_qstat)) 
  
  if len(result_qstat) == 0:
    qstat_out = 0
    logging.debug("check_qstat: Task "+str(id_in)+" is dead (unlisted).\n")
  elif 'finished' in result_qstat[0]:
    qstat_out = 0
    logging.debug("check_qstat: Task "+str(id_in)+" is dead (C).\n")
  elif ' C ' in result_qstat[2]:
    qstat_out = 0
    logging.debug("check_qstat: Task "+str(id_in)+" is dead (C).\n")
  elif ' Q ' in result_qstat[2]:
    qstat_out = 2
    logging.debug("check_qstat: Task "+str(id_in)+" is waiting (Q).")
  elif ' R ' in result_qstat[2]:
    qstat_out = 1
    logging.debug("check_qstat: Task "+str(id_in)+" is running (R).")
  else:
    qstat_out = 1
    logging.debug("check_qstat: Task "+str(id_in)+" is alive (unclear).")

  return qstat_out

#pbs_stat = check_qstat(110113)
#print("pbs_stat = "+str(pbs_stat))



for icase in range(icase_s, icase_e + 1):
  dir_work = dir_parent + "/" + dir_pre + str(icase)
  os.system("cd "+dir_parent+"/"+dir_pre+str(icase)+"/")
  
  logging.info("Calling qsub")
  #pbs_id = os.system("cd "+dir_parent+"/"+dir_pre+str(icase)+"/ && ""qsub -V " + pbs_script)
  cmd_tmp = "cd "+dir_parent+"/"+dir_pre+str(icase)+"/ && qsub -N "+pbs_pre+str(icase) + pbs_email + " -V "+pbs_script
  logging.debug("Executing:\n"+cmd_tmp)
  pbs_id = os.popen(cmd_tmp).readlines()[0]
  pbs_id = pbs_id.rstrip()
  logging.info('pbs_id = '+str(pbs_id))
  t_sub_submit = datetime.datetime.now()
  
  i_qstat = 1 
  n_qstat_q = 0
  n_qstat_r = 0
  while i_qstat > 0:
    i_qstat = check_qstat(pbs_id) ## return 0 (C), 1 (R), 2(Q)

    if i_qstat == 2:
      n_qstat_q = n_qstat_q + 1
      logging.debug("Get qstat Q for {0:d} time ".format(n_qstat_q))
      time.sleep(sleep_later)
    elif i_qstat == 1:
      n_qstat_r = n_qstat_r + 1
      logging.debug("Get qstat R for {0:d} time ".format(n_qstat_r))
      if n_qstat_r == 1:
        t_sub_begin = datetime.datetime.now()
        time.sleep(sleep_first)
      else:
        time.sleep(sleep_later)

  logging.info("Last task is finished. "+str(icase-icase_s+1)+"/"+str(ncase)+": "+str(pbs_id))
  t_sub_end = datetime.datetime.now()
  t_sub_queue = t_sub_begin - t_sub_submit
  t_sub_elap = t_sub_end - t_sub_begin

  t_elap = t_sub_end - t_begin
  t_sub_elap_mean = t_elap / (icase - icase_s + 1)
  t_end_est = t_begin + t_sub_elap_mean * ncase

  logging.info("Last task used "+str(t_sub_elap)+" in R, and "+str(t_sub_queue)+" in Q")
  logging.info("Average elapsed time for a task is "+str(t_sub_elap_mean))
  logging.info("Update finish-time estimation: "+str(t_end_est))
  
  sleep_first_old = sleep_first
  sleep_first = t_sub_elap.seconds + 24*3600* t_sub_elap.days - sleep_later
  if sleep_first < 0:
    sleep_first = sleep_later
  logging.debug("Update sleep_first from "+str(sleep_first_old)+" to "+str(sleep_first)+" seconds")

  logging.info("Copying restart files")
  cmd_tmp = "cd "+dir_parent+"/"+dir_pre+str(icase)+"/ && " + "sh copy_to_template.sh"
  logging.debug("Executing:\n "+cmd_tmp)
  os.system(cmd_tmp)
  
  cmd_tmp = "mv case_template "+dir_parent+"/"+dir_pre+str(icase+1)
  logging.debug("Executing:\n"+cmd_tmp)
  os.system(cmd_tmp)
  

