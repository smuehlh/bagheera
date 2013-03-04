require 'spec_helper'

describe PredictionsController do

	describe "Prediction page" do

		subject { page }
		before {visit root_path}
		it {should have_selector('title', text: full_title('Prediction'))}
		it {should have_selector('h1', text: 'Prediction')}
		it {should have_selector('h2', text: 'Genome file upload')}
		it {should have_selector('input', :class => 'ajax_upload_field')}

		it "should have the right links on the layout" do
		    click_link "Help"
		    page.should have_selector 'title', text: full_title('Help')
		    click_link "Contact"
		    page.should have_selector 'title', text: full_title('Contact')
		    click_link "Team"
		    page.should have_selector 'title', text: full_title('Team')
		    click_link "Prediction"
		    page.should have_selector 'title', text: full_title('Prediction')
		end

		it "has no file uploaded" do
			page.should have_content('No file uploaded')
			page.should have_selector('input', :id => 'predict_button')
			page.should have_selector('h1', text: 'Results', :visible => false)
			page.should have_selector('h2', text: 'Gene prediction and alignment options', :visible => false)
		end

		it "should upload file" do
			post upload_file_path, :uploaded_file => Rack::Test::UploadedFile.new(Rails.root.join('spec', 'fixtures', 'files', 'Candida_albicans_WO_1.fasta'), 'text/plain')
			response.should be_success
			response.should render_template('upload_file_ajax')
		end	
	end


	describe "Prediction workflow" do
		subject { page }

		it "has file uploaded" do
			get root_path
			get load_example_path # use get to keep one session
			render_template('upload_file_ajax')
			response.body.should have_content('Candida_albicans_WO_1.fasta')
			response.body.should have_selector('input', :id => "predict_button")
			response.body.should have_content("Alignment method")
			post predict_genes_path, algo: "gotoh", config: "tttt", species: "candida_albicans"
			render_template('predict_genes')

		end

	end
end
