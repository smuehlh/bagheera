source 'https://rubygems.org'

gem 'rails', '3.2.12'

# Bundle edge Rails instead:
# gem 'rails', :git => 'git://github.com/rails/rails.git'

group :development, :test do
  gem 'rspec-rails', '2.11.0' # test evironment
  gem 'guard-rspec', '1.2.1' # listen changes
  gem 'guard-spork', '0.3.2'
  gem 'spork', '0.9.2' # preload rails environment 
  gem "yard"
  gem "redcarpet"   # for yard to process .md files
  gem "rails_best_practices"
  gem 'simplecov', :require => false
  gem 'simplecov-rcov', :require => false
  gem "ci_reporter", require: false

end

# Gems used only for assets and not required
# in production environments by default.
group :assets do
  gem 'sass-rails',   '3.2.5'
  gem 'coffee-rails', '3.2.2'
  # gem 'therubyracer', :platforms => :ruby
  gem 'libv8', '~> 3.11.8.11'
  gem 'therubyracer', '~> 0.11.1', :platforms => :ruby
  gem 'uglifier', '>= 1.0.3'
end

gem 'jquery-rails'
gem 'jquery-fileupload-rails'
# gem 'jquery-validation-rails'

group :test do
  gem 'capybara', '1.1.2'
  gem 'rb-fsevent', '0.9.1', :require => false
  gem 'growl', '1.0.3'
  gem 'rb-inotify', '~> 0.8.8'
  gem 'factory_girl_rails', '4.1.0'
end

# To use ActiveModel has_secure_password
# gem 'bcrypt-ruby', '~> 3.0.0'

# To use Jbuilder templates for JSON
# gem 'jbuilder'

# Use unicorn as the app server
# gem 'unicorn'

# Deploy with Capistrano
# gem 'capistrano'

gem 'thin'
gem 'debugger'
gem "dalli"
gem 'fancybox-rails'
gem "peach"