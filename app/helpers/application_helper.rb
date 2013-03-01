module ApplicationHelper

  # Returns the full title on a per-page basis
  # @param page_title [String] Title of page
  # @return [String] Full tile of page (including "Bagheera | ")
  def full_title(page_title)
    base_title = "Bagheera"
    if page_title.empty?
      base_title
    else
      "#{base_title} | #{page_title}"
    end
  end
end
