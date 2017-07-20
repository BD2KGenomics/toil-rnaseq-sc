#!/usr/bin/env bash
echo "Log of Docker file?"
echo "Working directory:"
pwd
echo "List files:"
ls
echo "Running rscript:"
Rscript /image/ktsc.r --args "$@"
echo "Done."
