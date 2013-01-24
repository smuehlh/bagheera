require 'spec_helper'

describe "Predictions" do

	describe "Prediction page" do

		subject { page }
		before {visit root_path}
		it {should have_selector('title', text: full_title('Prediction'))}
		it {should have_selector('h1', text: 'Prediction')}

		it "should have the right links on the layout" do
		    visit root_path
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

	describe "file upload" do
		before {visit root_path}
		describe "when no file uploaded" do
			it {should_not have_selector('h1', text: 'Results')}
			it {should_not have_button("Predict")}
		end

		def genome_file_upload 
			genome_file = "#{Rails.root}/spec/fixtures/files/candida.fasta"
			Rack::Test::UploadedFile.new(genome_file, "text/plain")
		end
		describe "when POST to #upload_file" do
			before(:each) do
				xhr :post, upload_file: genome_file
				response.should be_success
			end
			it {should have_button("Predict")}
			it {should_not have_selector("h1", text: "Results")}
			it {should have_content("candida.fasta")}
		end


		# test for rendering partial
# it { should render_template(:partial => '_partialname') }
# xhr for xhttprequest
	end
end
