#!/usr/bin/env bash
echo "HELLO"
ls -R /data/
Rscript /image/ktsc.r --args "$@" 2>&1
