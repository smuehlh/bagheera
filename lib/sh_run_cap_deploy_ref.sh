# !/bin/sh

cd /fab8/4scp

echo "Prepare reference data"
mkdir bagheera_refdata
cp -r /tmp/cug/* bagheera_refdata

# delete all directories starting with numbers
find . -maxdepth 1 -type d -regextype sed -regex ".*/[0-9]*" -exec rm -rf {} \;

echo "Backup old reference data on cymo"

echo "Copy reference data to cymo"

# nicht noetig, weil auch ssh-keygen fuer fab8 & cymo reicht!!!
# -> klas fragen, ob ok wenn ich das einrichte


## add line to /etc/sudoers
## ANPASSEN: fab8 muss ohne password abfrage kopieren duerfen!
## fab8   ALL=(ALL:ALL) NOPASSWD: /path/to/your/script "" # "" to prevent user from running command with arbitray arguments


# file sh_run_cap_deploy.sh: run capistrano task to deploy reference data for bagheera on cymo

#cd /fab8/server/db_scripts/deployment/bagheera
#ap deploy_refdata
#exit 0
