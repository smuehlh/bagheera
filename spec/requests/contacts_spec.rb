require 'spec_helper'

describe "Contact page" do

	subject { page }
	before {visit contact_path}
	it {should have_selector('title', text: full_title('Contact'))}
	it {should have_selector('h1', text: 'Contact')}
	it {should have_selector('h2', text: 'Group leader')}
	it {should have_content("Dr. Martin Kollmar")}
	it {should have_content("mako@nmr.mpibpc.mpg.de")}

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