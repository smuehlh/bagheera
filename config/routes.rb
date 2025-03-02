Bagheera::Application.routes.draw do
 if ENV && ENV["PWD"] && (ENV["PWD"].include?("fab8") || ENV["PWD"].include?("jenkins")) then
    root to: 'static_pages#home'

    match '/prediction', to: 'predictions#search', as: 'prediction', :via => [:get]
    match '/team', to: 'static_pages#team', as: 'team', :via => [:get]
    match '/help', to: 'static_pages#help', as: 'help', :via => [:get]
    match '/contact', to: 'static_pages#contact', as: 'contact', :via => [:get]
    match '/specieslist', to: 'static_pages#specieslist', as: 'specieslist', :via => [:get]
    match '/translation', to: 'translations#start', as: 'translation', :via => [:get]

    match 'upload_file', to: 'predictions#upload_file', as: 'post/upload_file', :via => [:post]
    match 'upload_file', to: 'predictions#upload_file', as: 'get/upload_file', :via => [:get]
    match 'predict_genes', to: 'predictions#predict_genes', as: 'predict_genes', :via => [:post]
    match 'predict_more', to: 'predictions#predict_more', as: 'predict_more', :via => [:get]
    match 'show_alignment', to: 'predictions#show_alignment', as: 'show_alignment', :via => [:get]
    match 'calc_tree', to: 'tree#calc_tree', as: 'calc_tree', :via => [:post]
    match 'download', to: 'tree#download', as: 'download', :via => [:get]
    match 'transl_mrna', to: 'translations#transl_mrna', as: 'transl_mrna', :via => [:post]
    match 'transl_protein', to: 'translations#transl_protein', as: 'transl_protein', :via => [:post]
    match 'download_seq', to: 'translations#download', as: 'download_seq', :via => [:get]
    match 'upload_example', to: 'translations#upload_example', as: 'upload_example', :via => [:get]

    # match 'read_status', to: 'predictions#read_status', as: 'read_status', :via => [:get]
  else
    scope '/bagheera' do

      root to: 'static_pages#home'

      match '/prediction', to: 'predictions#search', as: 'prediction', :via => [:get]
      match '/team', to: 'static_pages#team', as: 'team', :via => [:get]
      match '/help', to: 'static_pages#help', as: 'help', :via => [:get]
      match '/contact', to: 'static_pages#contact', as: 'contact', :via => [:get]
      match '/specieslist', to: 'static_pages#specieslist', as: 'specieslist', :via => [:get]
      match '/translation', to: 'translations#start', as: 'translation', :via => [:get]

      match 'upload_file', to: 'predictions#upload_file', as: 'post/upload_file', :via => [:post]
      match 'upload_file', to: 'predictions#upload_file', as: 'get/upload_file', :via => [:get]
      match 'predict_genes', to: 'predictions#predict_genes', as: 'predict_genes', :via => [:post]
      match 'predict_more', to: 'predictions#predict_more', as: 'predict_more', :via => [:get]
      match 'show_alignment', to: 'predictions#show_alignment', as: 'show_alignment', :via => [:get]
      match 'calc_tree', to: 'tree#calc_tree', as: 'calc_tree', :via => [:post]
      match 'download', to: 'tree#download', as: 'download', :via => [:get]
      match 'transl_mrna', to: 'translations#transl_mrna', as: 'transl_mrna', :via => [:post]
      match 'transl_protein', to: 'translations#transl_protein', as: 'transl_protein', :via => [:post]
      match 'download_seq', to: 'translations#download', as: 'download_seq', :via => [:get]
      match 'upload_example', to: 'translations#upload_example', as: 'upload_example', :via => [:get]
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
