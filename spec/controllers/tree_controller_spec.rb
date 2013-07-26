require 'spec_helper'

describe TreeController do

  describe "GET 'calc_tree'" do
    it "returns http success" do
      get 'calc_tree'
      response.should be_success
    end
  end

  describe "GET 'gblocks'" do
    it "returns http success" do
      get 'gblocks'
      response.should be_success
    end
  end

  describe "GET 'fasttree'" do
    it "returns http success" do
      get 'fasttree'
      response.should be_success
    end
  end

end
