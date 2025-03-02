require 'spec_helper'

describe "StaticPages" do
  subject { page }

  shared_examples_for "all static pages" do
    it { should have_selector('h1',    text: heading) }
    it { should have_selector('title', text: full_title(page_title)) }
  end

  describe "Help page" do
    before { visit help_path } 
    let(:heading) {'Help'}
    let(:page_title) {'Help'}
    it_should_behave_like "all static pages"
  end

  describe "Team page" do
    before { visit team_path }
    let(:heading) {'Team'}
    let(:page_title) {'Team'}
    it_should_behave_like "all static pages"
  end

  describe "Contact page" do
    before { visit contact_path }
    let(:heading) {'Contact'}
    let(:page_title) {'Contact'}
    it_should_behave_like "all static pages"
  end

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
    click_link "Bagheera Team #{Time.now.year}"
    page.should have_selector 'title', text: full_title('Team')
    # test for links only the presence, not if they work correctly
    page.should have_link 'Impressum'
    page.should have_link 'Disney'
    page.should have_link 'link to diark'
    page.should have_link 'link to cymobase'
    page.should have_link 'link to motorprotein'
    page.should have_link 'link to MPG'
    page.should have_link 'MPI for biophysical chemistry'
  end
end
