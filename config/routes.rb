Bagheera::Application.routes.draw do
 if ENV && ENV["PWD"] && (ENV["PWD"].include?("fab8") || ENV["PWD"].include?("jenkins")) then
    root to: 'static_pages#home'

    match '/prediction', to: 'predictions#search', as: 'prediction', :via => [:get]
    match '/team', to: 'static_pages#team', as: 'team', :via => [:get]
    match '/help', to: 'static_pages#help', as: 'help', :via => [:get]
    match '/contact', to: 'static_pages#contact', as: 'contact', :via => [:get]

    match 'upload_file', to: 'predictions#upload_file', as: 'upload_file', :via => [:post]
    match 'predict_genes', to: 'predictions#predict_genes', as: 'predict_genes', :via => [:post]
    match 'load_example', to: 'predictions#load_example', as: 'load_example', :via => [:get]
    match 'show_alignment', to: 'predictions#show_alignment', as: 'show_alignment', :via => [:get]
    match 'predict_more', to: 'predictions#predict_more', as: 'predict_more', :via => [:get]
    match 'calc_tree', to: 'tree#calc_tree', as: 'calc_tree', :via => [:post]
    match 'download', to: 'tree#download', as: 'download', :via => [:get]

    # match 'read_status', to: 'predictions#read_status', as: 'read_status', :via => [:get]
    # match 'eval', to: 'predictions#mean_protein_length', as: 'eval', :via => [:get]
  else
    scope '/bagheera' do

      root to: 'predictions#search'

      match '/team', to: 'static_pages#team', as: 'team', :via => [:get]
      match '/help', to: 'static_pages#help', as: 'help', :via => [:get]
      match '/contact', to: 'static_pages#contact', as: 'contact', :via => [:get]

      match 'upload_file', to: 'predictions#upload_file', as: 'upload_file', :via => [:post]
      match 'predict_genes', to: 'predictions#predict_genes', as: 'predict_genes', :via => [:post]
      match 'load_example', to: 'predictions#load_example', as: 'load_example', :via => [:get]
      match 'show_alignment', to: 'predictions#show_alignment', as: 'show_alignment', :via => [:get]
      match 'predict_more', to: 'predictions#predict_more', as: 'predict_more', :via => [:get]
      match 'calc_tree', to: 'tree#calc_tree', as: 'calc_tree', :via => [:post]
      match 'download', to: 'tree#download', as: 'download', :via => [:get]
      # match 'eval', to: 'predictions#stat_conserved_pos', as: 'eval', :via => [:get]
    end
  end



  # The priority is based upon order of creation:
  # first created -> highest priority.

  # Sample of regular route:
  #   match 'products/:id' => 'catalog#view'
  # Keep in mind you can assign values other than :controller and :action

  # Sample of named route:
  #   match 'products/:id/purchase' => 'catalog#purchase', :as => :purchase
  # This route can be invoked with purchase_url(:id => product.id)

  # Sample resource route (maps HTTP verbs to controller actions automatically):
  #   resources :products

  # Sample resource route with options:
  #   resources :products do
  #     member do
  #       get 'short'
  #       post 'toggle'
  #     end
  #
  #     collection do
  #       get 'sold'
  #     end
  #   end

  # Sample resource route with sub-resources:
  #   resources :products do
  #     resources :comments, :sales
  #     resource :seller
  #   end

  # Sample resource route with more complex sub-resources
  #   resources :products do
  #     resources :comments
  #     resources :sales do
  #       get 'recent', :on => :collection
  #     end
  #   end

  # Sample resource route within a namespace:
  #   namespace :admin do
  #     # Directs /admin/products/* to Admin::ProductsController
  #     # (app/controllers/admin/products_controller.rb)
  #     resources :products
  #   end

  # You can have the root of your site routed with "root"
  # just remember to delete public/index.html.
  # root :to => 'welcome#index'

  # See how all your routes lay out with "rake routes"

  # This is a legacy wild controller route that's not recommended for RESTful applications.
  # Note: This route will make all actions in every controller accessible via GET requests.
  # match ':controller(/:action(/:id))(.:format)'
end
