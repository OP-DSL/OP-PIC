Search.setIndex({"alltitles":{"(1) Initialisation and Termination":[[0,"initialisation-and-termination"]],"(2) Dataset Layout":[[0,"dataset-layout"]],"(3) Parallel Loops":[[0,"parallel-loops"]],"(4) Particle injections":[[0,"particle-injections"]],"(5) HDF5 I/O":[[0,"hdf5-i-o"]],"(a) Allocate space for particles":[[1,"a-allocate-space-for-particles"]],"(a) Declare sets":[[1,"a-declare-sets"]],"(a) Direct loop":[[1,"a-direct-loop"]],"(a) Distributing data over MPI ranks":[[1,"a-distributing-data-over-mpi-ranks"]],"(b) Declare maps":[[1,"b-declare-maps"]],"(b) Indirect loop (single indirection)":[[1,"b-indirect-loop-single-indirection"]],"(b) Initialize the injected particles":[[1,"b-initialize-the-injected-particles"]],"(b) Partitioning data over MPI ranks":[[1,"b-partitioning-data-over-mpi-ranks"]],"(c) Declare data":[[1,"c-declare-data"]],"(c) Double Indirect loop":[[1,"c-double-indirect-loop"]],"(d) Declare constants":[[1,"d-declare-constants"]],"Compile OP-PIC library":[[2,"compile-op-pic-library"]],"Compile application":[[2,"compile-application"]],"Contents:":[[3,null]],"Dependencies":[[2,"dependencies"]],"Developer Guide":[[1,null]],"Example Application":[[1,"example-application"]],"Getting Started":[[2,null]],"Introduction":[[4,null]],"Licencing":[[4,"licencing"]],"OP-PIC C++ API":[[0,null]],"Original - Load mesh and initialization":[[1,"original-load-mesh-and-initialization"]],"Overview":[[4,"overview"]],"Run the application":[[2,"run-the-application"]],"Setup translator python environment (One time process)":[[2,"setup-translator-python-environment-one-time-process"]],"Step 1 - Preparing to use OP-PIC":[[1,"step-1-preparing-to-use-op-pic"]],"Step 2 - OP-PIC Declaration":[[1,"step-2-op-pic-declaration"]],"Step 3 - Parallel loop : opp_par_loop":[[1,"step-3-parallel-loop-opp-par-loop"]],"Step 4 - Move loop : opp_particle_move":[[1,"step-4-move-loop-opp-particle-move"]],"Step 5 - Particle injections":[[1,"step-5-particle-injections"]],"Step 6 - Global reductions":[[1,"step-6-global-reductions"]],"Step 7 - Handing it all to OP-PIC":[[1,"step-7-handing-it-all-to-op-pic"]],"Step 8 - Code generation":[[1,"step-8-code-generation"]],"Welcome to OP-PIC documentation!":[[3,null]]},"docnames":["api","developer_guide","getting_started","index","introduction"],"envversion":{"sphinx":65,"sphinx.domains.c":3,"sphinx.domains.changeset":1,"sphinx.domains.citation":1,"sphinx.domains.cpp":9,"sphinx.domains.index":1,"sphinx.domains.javascript":3,"sphinx.domains.math":2,"sphinx.domains.python":4,"sphinx.domains.rst":2,"sphinx.domains.std":2},"filenames":["api.rst","developer_guide.rst","getting_started.rst","index.rst","introduction.rst"],"indexentries":{"opp_arg_dat (c function)":[[0,"c.opp_arg_dat",false]],"opp_arg_gbl (c function)":[[0,"c.opp_arg_gbl",false]],"opp_decl_const (c function)":[[0,"c.opp_decl_const",false]],"opp_decl_dat (c function)":[[0,"c.opp_decl_dat",false]],"opp_decl_dat_hdf5 (c function)":[[0,"c.opp_decl_dat_hdf5",false]],"opp_decl_map (c function)":[[0,"c.opp_decl_map",false]],"opp_decl_map_hdf5 (c function)":[[0,"c.opp_decl_map_hdf5",false]],"opp_decl_particle_set (c function)":[[0,"c.opp_decl_particle_set",false]],"opp_decl_particle_set_hdf5 (c function)":[[0,"c.opp_decl_particle_set_hdf5",false]],"opp_decl_set (c function)":[[0,"c.opp_decl_set",false]],"opp_decl_set_hdf5 (c function)":[[0,"c.opp_decl_set_hdf5",false]],"opp_exit (c function)":[[0,"c.opp_exit",false]],"opp_inc_part_count_with_distribution (c function)":[[0,"c.opp_inc_part_count_with_distribution",false]],"opp_increase_particle_count (c function)":[[0,"c.opp_increase_particle_count",false]],"opp_init (c function)":[[0,"c.opp_init",false]],"opp_par_loop (c function)":[[0,"c.opp_par_loop",false]],"opp_particle_move (c function)":[[0,"c.opp_particle_move",false]]},"objects":{"":[[0,0,1,"c.opp_arg_dat","opp_arg_dat"],[0,0,1,"c.opp_arg_gbl","opp_arg_gbl"],[0,0,1,"c.opp_decl_const","opp_decl_const"],[0,0,1,"c.opp_decl_dat","opp_decl_dat"],[0,0,1,"c.opp_decl_dat_hdf5","opp_decl_dat_hdf5"],[0,0,1,"c.opp_decl_map","opp_decl_map"],[0,0,1,"c.opp_decl_map_hdf5","opp_decl_map_hdf5"],[0,0,1,"c.opp_decl_particle_set","opp_decl_particle_set"],[0,0,1,"c.opp_decl_particle_set_hdf5","opp_decl_particle_set_hdf5"],[0,0,1,"c.opp_decl_set","opp_decl_set"],[0,0,1,"c.opp_decl_set_hdf5","opp_decl_set_hdf5"],[0,0,1,"c.opp_exit","opp_exit"],[0,0,1,"c.opp_inc_part_count_with_distribution","opp_inc_part_count_with_distribution"],[0,0,1,"c.opp_increase_particle_count","opp_increase_particle_count"],[0,0,1,"c.opp_init","opp_init"],[0,0,1,"c.opp_par_loop","opp_par_loop"],[0,0,1,"c.opp_particle_move","opp_particle_move"]],"opp_arg_dat":[[0,1,1,"c.opp_arg_dat","acc"],[0,1,1,"c.opp_arg_dat","dat"]],"opp_arg_gbl":[[0,1,1,"c.opp_arg_gbl","acc"],[0,1,1,"c.opp_arg_gbl","data"],[0,1,1,"c.opp_arg_gbl","dim"],[0,1,1,"c.opp_arg_gbl","type"]],"opp_decl_const":[[0,1,1,"c.opp_decl_const","data"],[0,1,1,"c.opp_decl_const","dim"],[0,1,1,"c.opp_decl_const","name"]],"opp_decl_dat":[[0,1,1,"c.opp_decl_dat","data"],[0,1,1,"c.opp_decl_dat","dim"],[0,1,1,"c.opp_decl_dat","dtype"],[0,1,1,"c.opp_decl_dat","name"],[0,1,1,"c.opp_decl_dat","set"]],"opp_decl_dat_hdf5":[[0,1,1,"c.opp_decl_dat_hdf5","dim"],[0,1,1,"c.opp_decl_dat_hdf5","dtype"],[0,1,1,"c.opp_decl_dat_hdf5","file"],[0,1,1,"c.opp_decl_dat_hdf5","name"],[0,1,1,"c.opp_decl_dat_hdf5","set"]],"opp_decl_map":[[0,1,1,"c.opp_decl_map","dim"],[0,1,1,"c.opp_decl_map","from"],[0,1,1,"c.opp_decl_map","imap"],[0,1,1,"c.opp_decl_map","name"],[0,1,1,"c.opp_decl_map","to"]],"opp_decl_map_hdf5":[[0,1,1,"c.opp_decl_map_hdf5","dim"],[0,1,1,"c.opp_decl_map_hdf5","file"],[0,1,1,"c.opp_decl_map_hdf5","from"],[0,1,1,"c.opp_decl_map_hdf5","name"],[0,1,1,"c.opp_decl_map_hdf5","to"]],"opp_decl_particle_set":[[0,1,1,"c.opp_decl_particle_set","cells_set"],[0,1,1,"c.opp_decl_particle_set","name"]],"opp_decl_particle_set_hdf5":[[0,1,1,"c.opp_decl_particle_set_hdf5","cells_set"],[0,1,1,"c.opp_decl_particle_set_hdf5","file"],[0,1,1,"c.opp_decl_particle_set_hdf5","name"]],"opp_decl_set":[[0,1,1,"c.opp_decl_set","name"],[0,1,1,"c.opp_decl_set","size"]],"opp_decl_set_hdf5":[[0,1,1,"c.opp_decl_set_hdf5","file"],[0,1,1,"c.opp_decl_set_hdf5","name"]],"opp_inc_part_count_with_distribution":[[0,1,1,"c.opp_inc_part_count_with_distribution","p_set"],[0,1,1,"c.opp_inc_part_count_with_distribution","part_dist"],[0,1,1,"c.opp_inc_part_count_with_distribution","parts_to_insert"]],"opp_increase_particle_count":[[0,1,1,"c.opp_increase_particle_count","p_set"],[0,1,1,"c.opp_increase_particle_count","parts_to_insert"]],"opp_init":[[0,1,1,"c.opp_init","argc"],[0,1,1,"c.opp_init","argv"]],"opp_par_loop":[[0,1,1,"c.opp_par_loop","iter_type"],[0,1,1,"c.opp_par_loop","kernel"],[0,1,1,"c.opp_par_loop","name"],[0,1,1,"c.opp_par_loop","set"]],"opp_particle_move":[[0,1,1,"c.opp_particle_move","c2c_map"],[0,1,1,"c.opp_particle_move","kernel"],[0,1,1,"c.opp_particle_move","name"],[0,1,1,"c.opp_particle_move","p2c_map"],[0,1,1,"c.opp_particle_move","set"]]},"objnames":{"0":["c","function","C function"],"1":["c","functionParam","C function parameter"]},"objtypes":{"0":"c:function","1":"c:functionParam"},"terms":{"":1,"0":[0,1,2],"000":1,"03":1,"0th":1,"1":[2,3],"10":[1,2],"11":1,"12":1,"16":1,"1d":[0,1],"1st":1,"2":[2,3],"2057699":1,"21":1,"2d":[0,1],"2nd":1,"3":[2,3,4],"31":1,"32":2,"3d":[0,1],"3rd":1,"4":[2,3],"40":1,"47":1,"4th":1,"5":[2,3],"50":1,"54":1,"5th":1,"6":[2,3],"64":1,"7":3,"76":1,"8":[2,3],"85":1,"89":2,"9":1,"95":1,"A":[0,1,2],"As":1,"At":[0,1],"But":1,"By":0,"For":[0,1,2],"If":[0,1],"In":[0,1,2],"It":1,"No":1,"ONE":1,"One":[1,3],"That":1,"The":[0,1,2,4],"Then":[0,1],"There":1,"These":[1,4],"To":[0,1,2],"_cg":2,"ab":1,"abl":2,"about":1,"abov":1,"acc":0,"access":[0,1],"accord":1,"action":1,"ad":0,"add":1,"addit":[0,1,2],"addition":[0,1,2],"additon":1,"address":1,"adjac":0,"after":[0,1],"again":1,"algorithm":[1,4],"all":[0,3],"allow":1,"allthough":1,"almost":1,"along":[0,1],"alreadi":[0,1,2],"also":[0,1],"altern":0,"although":2,"alwai":1,"amd":2,"amount":0,"an":[0,1,2,4],"ani":[1,2],"anoth":1,"anticip":1,"ao":0,"api":[1,2,3,4],"app_":2,"app_cabanap":2,"app_femp":1,"app_fempic_opphc":1,"app_handcod":1,"app_nam":2,"applic":[0,3,4],"application_cpp_fil":[1,2],"approach":1,"appropri":[1,2],"ar":[0,1,2,4],"architectur":[2,4],"arg":1,"argc":[0,1],"argument":[0,1],"argv":[0,1],"ariti":1,"arrai":[0,1],"assign":0,"associ":0,"assum":1,"attach":[0,1],"automat":[2,4],"avail":[0,1],"avoid":1,"axi":0,"back":[0,1],"backend":2,"base":1,"basic":1,"batch":2,"becom":0,"been":[0,1],"befor":[0,1],"begin":1,"being":[0,2],"below":[0,1],"benefici":1,"benefit":1,"better":1,"between":[0,1],"bin":2,"bit":2,"blob":1,"block":1,"bodi":1,"bookkeep":1,"bool":1,"both":1,"bound":[0,1],"boundari":[0,1],"bounding_box":[0,1],"boundingbox":[0,1],"box":[0,1],"bsd":4,"buffer":0,"build":2,"c":[2,3,4],"c2c":1,"c2c_map":[0,1],"c2n_map":1,"c_centroid":1,"c_col":1,"c_color":1,"c_det":1,"c_ef":1,"c_gbl_id":[0,1],"c_id":1,"c_sd":1,"c_shape_deri":1,"c_to_c":1,"c_to_n":1,"c_v_c_map":1,"c_v_n_map":1,"c_vol":1,"c_volum":1,"cabana":2,"cabanap":2,"cach":1,"calc_pos_vel_kernel":1,"calcul":[0,1],"calculate_new_pos_vel":1,"call":[0,1],"can":[0,1,2],"capabl":1,"capac":0,"captur":1,"carri":1,"case":1,"cd":2,"cell":[0,1,2,4],"cell_det":1,"cell_ef":1,"cell_set":1,"cell_volum":1,"cells_set":0,"cgn":0,"chang":[1,2],"char":[0,1],"charg":1,"check":[0,1],"chosen":1,"clang":2,"class":1,"classifi":1,"claus":4,"cleanli":0,"close":1,"closer":1,"cluster":4,"code":[0,2,3,4],"coef":1,"coeff":1,"com":1,"command":[0,1,2],"comment":1,"commun":[0,1],"compat":1,"compil":[0,1,3],"complet":1,"comput":[0,1],"compute_ef_kernel":1,"compute_electric_field":1,"compute_ncd_kernel":1,"compute_node_charge_dens":1,"condit":1,"config":[1,2],"config_fil":2,"connect":1,"consid":1,"consist":[0,1],"const":[0,1],"const_":[0,1],"const_charg":1,"const_dt":1,"const_ion_veloc":1,"const_mass":1,"const_plasma_den":1,"const_spwt":1,"const_wall_potenti":1,"constant":0,"contain":[0,1,2],"contigu":0,"contribut":1,"control":0,"convert":1,"coordin":[0,1],"copi":[1,2],"core":[2,4],"correct":1,"could":1,"count":1,"cpp":[1,2],"cpu":[0,4],"crai":2,"creat":[0,1,2],"cu":1,"cuda":[1,2,4],"cuda_mpi":2,"current":[0,1,2,4],"custom":1,"dat":[0,1],"data":0,"datapoint":1,"dataset":3,"datatyp":0,"de":0,"decid":[0,1],"declar":[0,3],"default":0,"defin":[0,1],"definit":1,"degrad":1,"degre":1,"deletevalu":1,"demonstr":[1,2],"demostr":1,"densiti":1,"dep_charge_kernel":1,"depend":3,"deposit":1,"deposit_charge_on_nod":1,"deriv":1,"desc":1,"design":1,"destin":[0,1],"detail":[1,2,4],"develop":[0,3,4],"devic":1,"dh":1,"diagnost":0,"didxsize32":2,"differ":1,"dim":[0,1],"dimens":[0,1],"direct":0,"direct_hop":1,"directli":[0,1,2],"directori":[1,2],"disabl":1,"discuss":1,"distribut":[0,2,4],"distribute_data_over_rank":1,"do":1,"document":[1,4],"dof":1,"domain":[0,1,2,4],"done":1,"doubl":0,"download":2,"dp_rand":1,"dscotch_pthread":2,"dt":1,"dt_int":[0,1],"dt_real":[0,1],"dtype":0,"duct":1,"due":1,"dummi":1,"dummy_part_rand":1,"dummy_part_set":1,"dure":[0,1,2],"e":[1,2],"each":[0,1],"easi":0,"edsl":4,"effici":0,"either":0,"elabor":1,"electo":1,"electr":1,"electro":1,"electrostat":1,"element":[0,1],"elment":1,"els":1,"embed":[2,4],"enabl":[1,2],"end":[0,1],"endif":1,"enrich":1,"entail":1,"enum":0,"environ":[1,3],"equival":0,"error":0,"especi":1,"even":[0,1],"eventhough":1,"exampl":[0,2,3,4],"excalibur":1,"except":0,"execut":[0,1],"exist":0,"explain":1,"explan":2,"explicit":1,"explicitli":1,"extern":[0,1],"extract":[0,1],"face":1,"fact":1,"facto":0,"fals":1,"farawai":1,"fast":1,"fem":1,"fempic":1,"fempic_color_block":1,"fempic_convert_hdf5":1,"fempic_hdf5":1,"fempic_hdf5_opp":1,"fempic_misc_mesh_colour":1,"fempic_misc_mesh_load":1,"fempic_opp":1,"few":2,"field":1,"field_solv":1,"file":[0,1,2,4],"file_nam":1,"file_path":[1,2],"fill":1,"final":[0,1],"finalis":1,"find":[0,1],"finit":1,"first":1,"fix":1,"focu":1,"folder":[1,2],"follow":[1,2],"footprint":1,"form":[0,1],"format":[0,2],"found":[0,1,2],"fourth":1,"free":1,"freedom":1,"frequent":1,"from":[0,1,2],"full":1,"function":[0,1],"further":1,"g":[1,2],"gcc":2,"gener":[0,2,3,4],"geom":0,"geometr":0,"geomkwai":0,"get":[0,1,3],"get_final_max_valu":1,"get_final_max_values_kernel":1,"github":1,"give":1,"given":[0,1],"global":[0,3],"gnu":2,"go":1,"good":1,"gpu":[0,1,4],"grid":1,"grid_spac":[0,1],"guid":[0,3,4],"h":1,"ha":[0,1,2],"halo":[0,1],"hand":[2,3],"handcod":2,"handl":1,"have":[0,1,2],"hdf5":[1,2,3],"hdf5_install_path":2,"header":1,"henc":[0,1],"here":1,"high":[2,4],"higher":1,"hip":[1,2,4],"hip_mpi":2,"hold":1,"hole":1,"hop":[0,1],"how":[0,1],"howev":[0,1],"http":1,"i":[1,2,3,4],"id":0,"idea":1,"identifi":1,"idx":0,"if2c_map":1,"if2n_map":1,"if_area":1,"if_dist":1,"if_distrib":1,"if_n_po":1,"if_norm":1,"if_to_c":1,"if_to_n":1,"if_u_norm":1,"if_v_c_map":1,"if_v_n_map":1,"if_v_norm":1,"iface_area":1,"iface_dist":1,"iface_n_po":1,"iface_norm":1,"iface_set":1,"iface_u_norm":1,"iface_v_norm":1,"ifdef":1,"illustr":1,"imap":0,"immedi":1,"implement":[1,2],"implic":1,"impos":1,"includ":[0,1,2,4],"incorpor":1,"increas":[0,1],"increment":[0,1],"inde":1,"independ":1,"index":[0,1],"indic":[0,1],"indici":2,"indirect":0,"indirectli":1,"influenc":1,"inform":[0,1,4],"init_boundary_pot":1,"initi":[0,3],"initialis":[1,3],"initialz":1,"inject":3,"inject_ion":1,"inject_ions_kernel":1,"inlet":1,"inlet_faces_cel":1,"inlin":1,"inner":1,"input":0,"insert":[0,1],"instal":4,"instanc":1,"instead":0,"int":[0,1],"int_max":1,"integ":[0,1],"intel":2,"intend":1,"intern":1,"intervent":1,"introduc":1,"introduct":3,"invalid":1,"invoc":0,"invok":1,"ion":1,"ion_veloc":1,"issu":1,"iter":[0,1],"iter_typ":0,"its":[0,1],"kei":1,"kernel":[0,1],"keyword":0,"kind":1,"known":1,"kwai":0,"languag":[2,4],"larg":1,"last":[0,1],"later":[0,1],"layer":0,"layout":3,"leav":1,"length":1,"level":[2,4],"lib_nam":0,"libclang":2,"librari":[0,1,3],"licenc":3,"licens":4,"like":[0,1,2],"likewis":1,"limit":1,"line":[0,1],"linear":[1,2],"link":1,"list":0,"liter":[0,1],"load":3,"load_mesh":1,"loader":1,"log":1,"long":1,"longer":1,"look":1,"loop":3,"m":1,"magnet":1,"mai":[0,1,2],"main":[0,1],"mainli":[0,1],"maintain":1,"major":0,"make":[0,1,2],"manag":2,"mani":[2,4],"manual":[0,1,2],"map":0,"map1idx":1,"map2idx":1,"map3idx":1,"map4idx":1,"mappi":0,"march":1,"mark":1,"mass":1,"match":0,"matrix":[1,2],"max":1,"max_n_charge_den":1,"max_n_pot":1,"maxcoordin":[0,1],"maximum":[0,1],"mean":1,"mechan":1,"memori":[0,1,2,4],"mention":1,"mesh":[0,2,3,4],"mesh_cel":1,"mesh_nod":1,"method":1,"mh":1,"min_i":1,"min_lc":1,"mincoordin":[0,1],"mini":1,"minim":[0,1],"minimum":[0,1],"mix":1,"mode":1,"modifi":1,"more":[1,4],"most":1,"move":[0,3],"move_kernel":1,"move_particl":1,"movement":[0,1],"mover":1,"mpi":[0,2,4],"mpi_fin":0,"mpi_init":0,"mpicc":2,"mpicxx":2,"mpirun":2,"multi":[1,2,4],"multipl":1,"multipli":1,"must":0,"n_bnd_pot":1,"n_cell":1,"n_charge_den":1,"n_ifac":1,"n_ion_den":1,"n_node":1,"n_po":1,"n_pot":1,"n_pot0":1,"n_pot1":1,"n_pot2":1,"n_pot3":1,"n_potenti":1,"n_type":1,"n_vol":1,"n_volum":1,"name":[0,1,2],"ncd":1,"ncell":1,"need":[1,2],"neg":1,"neighbor":1,"neighbour":1,"neighour":0,"neptun":1,"new":1,"newli":1,"next":1,"nnode":1,"node":[0,1],"node_charge_den":1,"node_charge_den0":1,"node_charge_den1":1,"node_charge_den2":1,"node_charge_den3":1,"node_set":1,"node_volum":1,"non":0,"nonlinear":2,"note":1,"now":1,"nparticl":1,"nullptr":[0,1],"number":[0,1],"numer":1,"nv":1,"nvhpc":2,"nvidia":2,"o":[1,2,3],"object":0,"observ":1,"obtain":[1,2],"occur":[0,1],"offload":4,"omp":[1,2],"onc":1,"one":1,"onli":[0,1],"op":4,"op_read":1,"open":[1,4],"openmp":[2,4],"oper":1,"opp":[0,1],"opp_access":0,"opp_allocation_multipl":1,"opp_arg":0,"opp_arg_dat":[0,1],"opp_arg_gbl":[0,1],"opp_bool":1,"opp_c2c":1,"opp_dat":[0,1],"opp_data_typ":0,"opp_decl_const":[0,1],"opp_decl_dat":[0,1],"opp_decl_dat_hdf5":[0,1],"opp_decl_map":[0,1],"opp_decl_map_hdf5":[0,1],"opp_decl_particle_set":[0,1],"opp_decl_particle_set_hdf5":[0,1],"opp_decl_set":[0,1],"opp_decl_set_hdf5":[0,1],"opp_do_onc":1,"opp_exit":[0,1],"opp_inc":[0,1],"opp_inc_part_count_with_distribut":[0,1],"opp_increase_particle_count":[0,1],"opp_init":[0,1],"opp_init_direct_hop":[0,1],"opp_int":1,"opp_iterate_al":[0,1],"opp_iterate_inject":[0,1],"opp_iterate_typ":0,"opp_kernel":1,"opp_lib":2,"opp_map":[0,1],"opp_max":[0,1],"opp_min":0,"opp_move_particl":1,"opp_p2c":1,"opp_par_loop":[0,3],"opp_particle_done_mov":1,"opp_particle_mov":[0,3],"opp_particle_move_don":1,"opp_particle_need_mov":1,"opp_particle_need_remov":1,"opp_partit":[0,1],"opp_path":[1,2],"opp_point":[0,1],"opp_printf":1,"opp_profil":1,"opp_read":[0,1],"opp_real":1,"opp_root":1,"opp_run_on_root":1,"opp_rw":[0,1],"opp_set":[0,1],"opp_templ":1,"opp_transl":[1,2],"opp_venv":2,"opp_writ":[0,1],"optim":1,"option":[0,1,2],"order":1,"origin":3,"other":[0,1,2],"otherwis":0,"out":1,"outer":1,"outlin":1,"output":[0,1],"outsid":1,"over":0,"overal":1,"overhead":1,"overview":3,"own":1,"p2c":1,"p2c_map":[0,1],"p_lc":1,"p_po":1,"p_posit":1,"p_set":[0,1],"p_vel":1,"p_veloc":1,"packag":2,"page":[1,4],"parallel":[2,3,4],"parallelis":0,"param":2,"paramet":0,"parellelis":4,"parmeti":[0,2],"parmetis_geom":[0,1],"parmetis_install_path":2,"parmetis_kwai":[0,1],"part_dist":[0,1],"part_lc":1,"part_mesh_rel":1,"part_po":1,"part_posit":1,"part_vel":1,"part_veloc":1,"particl":[2,3,4],"particle_set":1,"particles_set":1,"particularli":[0,4],"partit":[0,2],"partitin":0,"partition":1,"parts_to_insert":[0,1],"pass":[0,1],"path":1,"pdf":1,"per":[0,1],"perform":1,"petsc":[0,1,2],"petsc_install_path":2,"petscfin":0,"petsciniti":0,"pic":4,"place":1,"plasma_den":1,"platform":2,"point":1,"point_lc":1,"point_po":1,"pointer":[0,1],"pos_dat":[0,1],"posit":[0,1],"possibl":0,"possibli":1,"potenti":1,"pre":1,"prepar":3,"present":[0,1],"previous":0,"primari":0,"prime_map":0,"prime_set":0,"print":1,"printf":1,"prior":[0,1,2],"probabl":1,"process":3,"processor":1,"program":[1,2],"project":4,"provid":[0,1,2,4],"python":3,"python3":[1,2],"random":1,"rank":0,"rate":1,"rather":1,"reach":1,"reachabl":0,"read":[0,1],"readi":2,"readm":[1,2],"realloc":1,"reduc":1,"reduct":[0,3],"regular":1,"relat":[0,1],"releas":4,"relev":1,"remain":1,"remov":1,"renumb":0,"report":1,"repositori":2,"request":1,"requir":[0,1,2],"resiz":1,"resolut":1,"retain":1,"return":[0,1],"rocm":2,"root":1,"routin":[0,1],"run":[0,1,3],"runtim":[0,1],"same":1,"save":1,"scale":1,"scatter":0,"scheme":[0,1],"scope":[0,1],"script":2,"search":1,"search_next_cel":1,"second":1,"section":[0,1],"see":[0,1,4],"seen":1,"select":1,"separ":1,"seq":[1,2],"sequenti":[1,2,4],"set":0,"set_capac":1,"setup":[1,3],"setup_venv":2,"seven":1,"sh":2,"shared_ptr":1,"shift":1,"should":[0,1],"show":1,"showcas":1,"shuffl":1,"signal":1,"significantli":1,"similar":1,"similarli":1,"simpli":1,"simul":[0,1,2],"sinc":[0,1],"singl":0,"six":1,"size":[0,1],"slurm":2,"so":1,"soa":0,"solv":1,"solver":[1,2],"some":[0,1,2],"sort":1,"sourc":[0,2,4],"space":0,"spars":[1,2],"special":1,"specif":[0,1,2,4],"specifi":[0,1,2],"spent":1,"spwt":1,"srun":2,"stage":1,"standard":0,"start":[1,3],"statement":1,"static":1,"statu":1,"std":[0,1],"step":3,"storag":1,"store":[0,1],"strategi":[0,1],"stream":1,"string":[0,1],"struct":0,"structur":[0,1,2],"subroutin":1,"suggest":0,"sum":0,"summar":1,"suppli":0,"support":[2,4],"surround":1,"switch":1,"synchron":1,"system":0,"t":0,"tabl":[0,1],"take":[0,1],"target":[0,1,4],"templat":1,"temporari":1,"termin":[1,3],"test":1,"tetrahedr":1,"text":1,"than":1,"thei":1,"them":1,"thi":[0,1,2],"third":1,"though":1,"thread":4,"three":1,"thrid":1,"through":[0,1],"thu":1,"time":[0,1,3],"timer":1,"tn":1,"togeth":1,"tool":2,"toolchain":2,"top":0,"toward":1,"track":1,"transform":0,"translat":[0,1,3],"true":1,"tutori":1,"two":1,"type":[0,1],"under":[0,1,4],"underli":0,"understand":1,"uniqu":1,"unless":0,"unnecessari":1,"unstructur":[1,2,4],"until":1,"us":[0,2,3,4],"usag":1,"use_mpi":1,"use_petsc":0,"user":[0,1],"util":1,"v":[1,2],"valid":[0,1],"valu":[0,1],"var_nam":0,"vari":1,"variabl":[0,1,2],"variant":4,"varieti":2,"variou":0,"veloc":1,"version":[1,2],"via":[1,2],"view":1,"void":[0,1],"volum":1,"wai":1,"walkthrough":4,"wall":1,"wall_potenti":1,"we":[0,1],"weight":1,"well":0,"were":1,"what":[0,1],"when":[0,1],"where":[0,1],"whether":0,"which":[0,1],"while":1,"wish":1,"wit":1,"within":[0,1],"without":[1,2],"word":1,"work":1,"would":0,"wrapper":2,"write":[0,1,2,4],"written":[1,2],"you":[1,2],"your":[1,2],"zero":[0,1]},"titles":["OP-PIC C++ API","Developer Guide","Getting Started","Welcome to OP-PIC documentation!","Introduction"],"titleterms":{"1":[0,1],"2":[0,1],"3":[0,1],"4":[0,1],"5":[0,1],"6":1,"7":1,"8":1,"One":2,"all":1,"alloc":1,"api":0,"applic":[1,2],"b":1,"c":[0,1],"code":1,"compil":2,"constant":1,"content":3,"d":1,"data":1,"dataset":0,"declar":1,"depend":2,"develop":1,"direct":1,"distribut":1,"document":3,"doubl":1,"environ":2,"exampl":1,"gener":1,"get":2,"global":1,"guid":1,"hand":1,"hdf5":0,"i":0,"indirect":1,"initi":1,"initialis":0,"inject":[0,1],"introduct":4,"layout":0,"librari":2,"licenc":4,"load":1,"loop":[0,1],"map":1,"mesh":1,"move":1,"mpi":1,"o":0,"op":[0,1,2,3],"opp_par_loop":1,"opp_particle_mov":1,"origin":1,"over":1,"overview":4,"parallel":[0,1],"particl":[0,1],"partit":1,"pic":[0,1,2,3],"prepar":1,"process":2,"python":2,"rank":1,"reduct":1,"run":2,"set":1,"setup":2,"singl":1,"space":1,"start":2,"step":1,"termin":0,"time":2,"translat":2,"us":1,"welcom":3}})