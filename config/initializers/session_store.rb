# Be sure to restart your server when you modify this file.

# Bagheera::Application.config.session_store :cookie_store, key: '_bagheera_session'
# Bagheera::Application.config.session_store :mem_cache_store, {:namespace => "Bagheera", :compression => true}

# Use the database for sessions instead of the cookie-based default,
# which shouldn't be used to store highly confidential information
# (create the session table with "rails generate session_migration")
# Bagheera::Application.config.session_store :active_record_store

require 'action_dispatch/middleware/session/dalli_store'
Bagheera::Application.config.session_store :dalli_store, memcache_server: ['localhost'], namespace: 'Bagheera', expire_after: 1.day, compress: true, value_max_bytes: 2.megabyte
