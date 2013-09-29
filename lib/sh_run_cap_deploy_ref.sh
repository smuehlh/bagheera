#!/bin/bash

# file sh_run_cap_deploy.sh: run capistrano task to deploy reference data for bagheera on cymo

cd /fab8/server/db_scripts/deployment/bagheera
cap deploy_refdata
exit 0
