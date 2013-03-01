# encoding: utf-8
require 'spec_helper'

describe "Team page" do

	subject { page }
	before {visit team_path}
	it {should have_selector('title', text: full_title('Team'))}
	it {should have_selector('h1', text: 'Team')}
	it {should have_selector('h2', text: 'Current Team Members')}
	it {should have_content("Dr. Martin Kollmar")}
	it {should have_content('Stefanie Mühlhausen')}
	it {should have_link ('group homepage')}

	it "should have the right links on the layout" do
	    visit contact_path
	    click_link "Help"
	    page.should have_selector 'title', text: full_title('Help')
	    click_link "Contact"
	    page.should have_selector 'title', text: full_title('Contact')
	    click_link "Team"
	    page.should have_selector 'title', text: full_title('Team')
	    click_link "Prediction"
	    page.should have_selector 'title', text: full_title('Prediction')
	end
end