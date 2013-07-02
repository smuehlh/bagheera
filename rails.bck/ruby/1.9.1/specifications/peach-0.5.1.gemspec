# -*- encoding: utf-8 -*-

Gem::Specification.new do |s|
  s.name = "peach"
  s.version = "0.5.1"

  s.required_rubygems_version = Gem::Requirement.new(">= 0") if s.respond_to? :required_rubygems_version=
  s.authors = ["Ben Hughes"]
  s.date = "2012-07-03"
  s.email = "ben@pixelmachine.org"
  s.homepage = "http://peach.rubyforge.org"
  s.require_paths = ["lib"]
  s.rubygems_version = "1.8.24"
  s.summary = "Parallel Each and other parallel things"

  if s.respond_to? :specification_version then
    s.specification_version = 3

    if Gem::Version.new(Gem::VERSION) >= Gem::Version.new('1.2.0') then
      s.add_development_dependency(%q<rake>, [">= 0"])
      s.add_development_dependency(%q<shoulda>, [">= 0"])
    else
      s.add_dependency(%q<rake>, [">= 0"])
      s.add_dependency(%q<shoulda>, [">= 0"])
    end
  else
    s.add_dependency(%q<rake>, [">= 0"])
    s.add_dependency(%q<shoulda>, [">= 0"])
  end
end
