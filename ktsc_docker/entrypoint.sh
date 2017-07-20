#!/usr/bin/env bash
pwd
ls
echo "Here are the args:$@"
Rscript /image/ktsc.r --args "$@" 2>&1
