require 'spec_helper'

describe "Help page" do

	subject { page }
	before {visit help_path}
	it {should have_selector('title', text: full_title('Help'))}
	it {should have_selector('h1', text: 'Help')}
	it {should have_content("file upload")}

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