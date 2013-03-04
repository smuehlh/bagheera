# -*- encoding: utf-8 -*-

Gem::Specification.new do |s|
  s.name = "code_analyzer"
  s.version = "0.3.1"

  s.required_rubygems_version = Gem::Requirement.new(">= 0") if s.respond_to? :required_rubygems_version=
  s.authors = ["Richard Huang"]
  s.date = "2012-12-18"
  s.description = "a code analyzer tool which extracted from rails_best_practices, it helps you easily build your own code analyzer tool."
  s.email = ["flyerhzm@gmail.com"]
  s.homepage = "https://github.com/flyerhzm/code_analyzer"
  s.require_paths = ["lib"]
  s.rubygems_version = "1.8.24"
  s.summary = "a code analyzer helps you build your own code analyzer tool."

  if s.respond_to? :specification_version then
    s.specification_version = 3

    if Gem::Version.new(Gem::VERSION) >= Gem::Version.new('1.2.0') then
      s.add_runtime_dependency(%q<sexp_processor>, [">= 0"])
      s.add_development_dependency(%q<rspec>, [">= 0"])
      s.add_development_dependency(%q<rake>, [">= 0"])
    else
      s.add_dependency(%q<sexp_processor>, [">= 0"])
      s.add_dependency(%q<rspec>, [">= 0"])
      s.add_dependency(%q<rake>, [">= 0"])
    end
  else
    s.add_dependency(%q<sexp_processor>, [">= 0"])
    s.add_dependency(%q<rspec>, [">= 0"])
    s.add_dependency(%q<rake>, [">= 0"])
  end
end
