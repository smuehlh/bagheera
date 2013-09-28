class ApplicationController < ActionController::Base
	protect_from_forgery

	# rescue_from ActionView::Template::Error do |exception|
	# 	render "error", formats: [:js]
	# end
end
