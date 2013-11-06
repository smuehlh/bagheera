# !/bin/sh

# copies Bagheera's reference data from fab8 to cymo
# do not run this script as root (as it will do ssh and scp), but as fab8
# fab8 has a ssh-key shared with wwwrun@cymo

# initialize variables
if ! [ $# -eq 1 ] ; then
	FAB8DATA='/fab8/server/bagheera/db'
else
	FAB8DATA=$1
fi
WEBAPPDATA='/usr/local/cymo/bagheera/db'
DIR='ref_data'
# reference me like ${DATA}

clean_up () {
	# housekeeping at program exit
	ssh wwwrun@cymo "cd ${WEBAPPDATA} && mv ${DIR}_old ${DIR}"
	exit 	
}

# save old reference data on cymo
ssh wwwrun@cymo "cd ${WEBAPPDATA} && mv ${DIR}/ ${DIR}_old"

# scp data to cymo
scp -rq ${FAB8DATA}/${DIR} wwwrun@cymo:${WEBAPPDATA}

# restore old reference data on cymo in case scp was not successful
STAT=$?
if [ $STAT -eq 0 ] ; then
	# housekeeping on success
	ssh wwwrun@cymo "cd ${WEBAPPDATA} && rm -rf ${DIR}_old"
	exit 0
else
	# no success 
	clean_up
fi

trap clean_up 1 2 3 13 15 