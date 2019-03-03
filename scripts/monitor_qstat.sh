#!/bin/bash
# Usage: bash monitor_qstat.sh <user>

N=5
USER=$1

while true; do
  clear
  # print current status
  qstat -u $USER
  
  echo -e "--------------\n"
  # print observed peak mem usage
  qstat -u $USER | grep -o "[0-9]*.torquelord" | cut -d'.' -f1 |
  xargs -I% sh -c \
      'echo Job % && echo $(qstat -f % | grep "resources_used.mem")'
  sleep $N
done
