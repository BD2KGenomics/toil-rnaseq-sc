#!/usr/bin/env bash
ls -R /data/
Rscript /image/ktsc.r --args "$@" 2>&1
