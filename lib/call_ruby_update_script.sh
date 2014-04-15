#!/bin/bash

# load rvm ruby
source /home/fab8/.rvm/environments/ruby-1.9.2-p290@bagheera
# call ruby script
ruby /fab8/server/bagheera/lib/update_ref_fab8.rb

sh /fab8/server/bagheera/lib/copy_bagheera_refdata_to_cymo.sh
