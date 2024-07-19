OP-PIC C++ API
==============

(1) Initialisation and Termination
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. c:function:: void opp_init(int argc, char **argv)

   This routine must be called before all other OP-PIC routines. 
   If `USE_PETSC` is defined during compilation, `opp_init` will initialize PETSc library as well (using :c:func:`PetscInitialize()`).
   Under MPI back-ends, this routine also calls :c:func:`MPI_Init()` unless its already called previously.

   :param argc: The number of command line arguments.
   :param argv: The command line arguments, as passed to :c:func:`main()`.

.. c:function:: void opp_exit()

   This routine must be called last to cleanly terminate the OP-PIC runtime. 
   If `USE_PETSC` is defined during compilation, `opp_exit` will finalize PETSc library as well (using :c:func:`PetscFinalize()`).
   Under MPI back-ends, this routine also calls :c:func:`MPI_Finalize()` unless its has been called previously. 
   A runtime error will occur if :c:func:`MPI_Finalize()` is called after :c:func:`opp_exit()`.

.. c:function:: opp_set opp_decl_set(int size, char const *name)

   This routine declares a mesh set.

   :param size: Number of set elements.
   :param name: A name to be used for output diagnostics.
   :returns: A set ID.

.. c:function:: opp_set opp_decl_particle_set(char const *name, opp_set cells_set)

   This routine declares a particle set with zero particles in it.

   :param name: A name to be used for output diagnostics.
   :param size: The underlying mesh set related.
   :returns: A set ID.

.. c:function:: opp_set opp_decl_particle_set(int size, char const *name, opp_set cells_set)

   This routine declares a particle set with a given number of particles in it.

   :param size: Number of particles in the set.
   :param name: A name to be used for output diagnostics.
   :param size: The underlying mesh set related.
   :returns: A set ID.

.. c:function:: opp_map opp_decl_map(opp_set from, opp_set to, int dim, int *imap, char const *name)

   This routine defines a mapping between sets. 
   
   :param from: Source set.
   :param to: Destination set.
   :param dim: Number of mappings per source element.
   :param imap: Mapping table.
   :param name: A name to be used for output diagnostics.

   .. note::
      imap should contain a valid data array pointer. The only exception is to provide `nullptr` given the from set is a particle set with zero particles in it.
   
.. c:function:: opp_dat opp_decl_dat(opp_set set, int dim, opp_data_type dtype, void *data, char const *name)

   This routine defines a dataset.

   :param set: The set the data is associated with (can be a mesh set or a particle set).
   :param dim: Number of data elements per set element.
   :param type: The datatype, This can be a type in :c:expr:`opp_data_type` enum (:c:expr:`DT_INT` or :c:expr:`DT_REAL`, other type may be added later).
   :param data: Input data of type :c:type:`T` (checked for consistency with **type** at run-time). The data must be provided in AoS form with each of the **dim** elements per set element contiguous in memory.
   :param name: A name to be used for output diagnostics.

   .. note::
      At present **dim** must be an integer literal or a :c:expr:`#define`

.. c:function:: void opp_decl_const(int dim, T* data, const char* name)

   This routine defines constant data with global scope that can be used in kernel functions.

   :param dim: Number of data elements. For maximum efficiency this should be an integer literal or a :c:expr:`#define`.
   :param data: A pointer to the data, checked for type consistency at run-time.
   :param name: The name as a string that the kernels access it, (should be :c:expr:`CONST_+<var_name>`).

   .. note::
      The variable is available in the kernel functions with type :c:expr:`T` with type :c:expr:`T*`. Hence even if **dim** is :c:expr:`1`, it should be accessed as :c:expr:`CONST_+<var_name>[0]` within the kernel.

.. c:function:: void opp_init_direct_hop(double grid_spacing, int dim, const opp_dat c_gbl_id, const opp::BoundingBox& bounding_box)

   Existance of this routine suggest the code-generator to extract information from :c:func:`opp_particle_move()` to generate direct-hop related code. Specifically, it will create direct-hop structure data and include direct hop related code within the particle move loop.

   :param grid_spacing: Direct hop structured mesh spacing.
   :param dim: Dimension of the simulation (1D, 2D or 3D).
   :param c_gbl_id: An `opp_dat` with global cell indices (mainly required to translate cell indices in an MPI code simulation, since mappings get renumbered).
   :param bounding_box: A `opp::BoundingBox` indicating the simulation boundaries

   .. note::
      The bounding box object can be created by providing a mesh dat that has its positions (like node positions) using :c:expr:`opp::BoundingBox(const opp_dat pos_dat, int dim)` or by providing the calculated minimum and maximum domain coordinates using :c:expr:`opp::BoundingBox(int dim, opp_point minCoordinate, opp_point maxCoordinate)`.

.. c:function:: void opp_partition(std::string lib_name, opp_set prime_set, opp_map prime_map, opp_dat dat)

   This routine controls the partitioning of the sets used for distributed memory parallel execution.

   :param lib_name: The partitioning library to use, see below.
   :param prime_set: Specifies the set to be partitioned.
   :param prime_map: Specifies the map to be used to create adjacency lists for the **prime_set**. Required if using :c:expr:`"KWAY"` or :c:expr:`"GEOMKWAY"` (defaulted to :c:expr:`nullptr`).
   :param dat: Specifies the :c:expr:`opp_dat` required for the partitioning strategy (defaulted to :c:expr:`nullptr`).

   The current options for **lib_name** are:

   - :c:expr:`"PARMETIS_KWAY"`: Uses the Kway routine of `ParMETIS <http://glaros.dtc.umn.edu/gkhome/metis/parmetis/overview>`_ library. Required to provide geometric coordinates of the **prime_set** through the data `opp_dat`.
   - :c:expr:`"PARMETIS_GEOM"`: Uses the Geom routine of `ParMETIS <http://glaros.dtc.umn.edu/gkhome/metis/parmetis/overview>`_ library. Required to provide a primary map.
   - :c:expr:`"EXTERNAL"`: External partitioning defining user specified partitioning scheme. The ranks for the primary set element to be, should be provided as a `opp_dat`. This routine can be used to partition the cells along the major particle movement axis to minimize particle communications.

   .. note::
      In addition to partitining the primary set, :c:func:`opp_partition` will partition the other mesh sets using `opp_maps` available. Then Halos and halo communication buffers will be created. Additionally, halo related particle communication data will get initialized.

(2) Dataset Layout
^^^^^^^^^^^^^^^^^^

By default OP-PIC stores data in CPUs as AoS (Array of Structs) layout, matching what is supplied to :c:func:`opp_decl_dat()` and :c:func:`opp_decl_map()`, however, on GPUs, OP-PIC use SoA (Struct of Arrays) layouts using transformations.

(3) Parallel Loops
^^^^^^^^^^^^^^^^^^

.. c:function:: void opp_par_loop(void (*kernel)(...), char const *name, opp_set set, opp_iterate_type iter_type, ...)

   This routine executes a parallelised loop over the given **set**, with arguments provided by the :c:func:`opp_arg_gbl()` and :c:func:`opp_arg_dat()` routines.
   When set is a particle set, it will make use of :c:expr:`iter_type` to decide whether to iterate over all particles or only over the injected particles.

   :param kernel: The kernel function to execute. The number of arguments to the kernel should match the number of :c:type:`opp_arg` arguments provided to this routine.
   :param name: A name to be used for output diagnostics.
   :param set: The set to loop over.
   :param iter_type: The iteration type. Possible values are, :c:expr:`OPP_ITERATE_ALL` and :c:expr:`OPP_ITERATE_INJECTED`
   :param ...: The :c:type:`opp_arg` arguments passed to each invocation of the kernel.

.. c:function:: void opp_particle_move(void (*kernel)(...), char const *name, opp_set set, opp_map c2c_map, opp_map p2c_map, ...)

   This routine executes a parallelised loop over the given **particle set**, with arguments provided by the :c:func:`opp_arg_gbl()` and :c:func:`opp_arg_dat()` routines.

   :param kernel: The kernel function to execute. The number of arguments to the kernel should match the number of :c:type:`opp_arg` arguments provided to this routine.
   :param name: A name to be used for output diagnostics.
   :param set: The set to loop over.
   :param c2c_map: A cell to cell mapping to find the neighour cell from the current cell.
   :param p2c_map: The particle to cell mapping to map to the underlying cell.
   :param ...: The :c:type:`opp_arg` arguments passed to each invocation of the kernel.

.. c:function:: opp_arg opp_arg_gbl(T *data, int dim, char *type, opp_access acc)

   This routine defines an :c:type:`opp_arg` that may be used either to pass non-constant read-only data or to compute a global sum, maximum or minimum.

   :param data: Source or destination data array.
   :param dim: Number of data elements.
   :param type: The datatype as a string. This is checked for consistency with **data** at run-time.
   :param acc: The access type.

   Valid access types for this routine are:

   - :c:data:`OPP_READ`: Read-only.
   - :c:data:`OPP_INC`: Global reduction to compute a sum.
   - :c:data:`OPP_MAX`: Global reduction to compute a maximum.
   - :c:data:`OPP_MIN`: Global reduction to compute a minimum.

**For direct arguments**

.. c:function:: opp_arg opp_arg_dat(opp_dat dat, opp_access acc) 

**For single indirect arguments**

.. c:function:: opp_arg opp_arg_dat(opp_dat dat, opp_map p2c_map, opp_access acc)
   opp_arg opp_arg_dat(opp_dat dat, int idx, opp_map map, opp_access acc)

**For double indirect arguments**

.. c:function:: opp_arg opp_arg_dat(opp_dat dat, int idx, opp_map map, opp_map p2c_map, opp_access acc) 

   This routine defines an :c:type:`opp_arg` that can be used to pass a dataset either directly attached to the target :c:type:`opp_set` or attached to an :c:type:`opp_set` reachable through a mapping.

   :param dat: The dataset.
   :param acc: The access type.
   :param p2c_map: Map from a particle to underlying cell.
   :param idx: The per-set-element index into the map to use. 
   :param map: The mapping to use for indirection.

   Valid access types for this routine are:

   - :c:data:`OPP_READ`: Read-only.
   - :c:data:`OPP_WRITE`: Write-only.
   - :c:data:`OPP_RW`: Read and write.
   - :c:data:`OPP_INC`: Increment or global reduction to compute a sum.

   .. note::
      An example of how these API calls are used in an application can be found in the Developer Guide Section.

(4) Particle injections
^^^^^^^^^^^^^^^^^^^^^^^

.. c:function:: void opp_increase_particle_count(opp_set p_set, int parts_to_insert)

   This routine will increase the particle set size and capacity by the given amount.

   :param p_set: Particle set.
   :param parts_to_insert: Number of particles to insert.

.. c:function:: void opp_inc_part_count_with_distribution(opp_set p_set, int parts_to_insert, opp_dat part_dist)

   This routine will increase the particle set size and capacity by the given amount and assigns :c:expr:`p2c_map` using :c:expr:`part_dist`.

   :param p_set: Particle set.
   :param parts_to_insert: Number of particles to insert.
   :param part_dist: Particle distribution :c:expr:`opp_dat`.

   .. note::
      An example of how these API calls are used in an application can be found in the Developer Guide Section.

(5) HDF5 I/O
^^^^^^^^^^^^

`HDF5 <https://www.hdfgroup.org/solutions/hdf5/>`_ has become the *de facto* format for parallel file I/O, with various other standards like `CGNS <https://cgns.github.io/hdf5.html>`_ layered on top. 
To make it as easy as possible for users to develop distributed-memory OP-PIC applications, we provide alternatives to some of the OP-PIC routines in which the data is read by OP-PIC from an HDF5 file, instead of being supplied by the user. 
This is particularly useful for distributed memory MPI systems where the user would otherwise have to manually scatter data arrays over nodes prior to initialisation.

.. c:function:: opp_set opp_decl_set_hdf5(char const *file, char const *name)

   Equivalent to :c:func:`opp_decl_set()` but takes a **file** instead of **size**, reading in the set size from the HDF5 file using the keyword **name**.

.. c:function:: opp_set opp_decl_particle_set_hdf5(char const *file, char const *name, opp_set cells_set)

   Equivalent to :c:func:`opp_decl_particle_set()` but takes a **file** instead of **size**, reading in the set size from the HDF5 file using the keyword **name**.

.. c:function:: opp_map opp_decl_map_hdf5(opp_set from, opp_set to, int dim, char const *file, char const *name)

   Equivalent to :c:func:`opp_decl_map()` but takes a **file** instead of **imap**, reading in the mappiing table from the HDF5 file using the keyword **name**.

.. c:function:: opp_dat opp_decl_dat_hdf5(opp_set set, int dim, opp_data_type dtype, char const *file, char const *name)

   Equivalent to :c:func:`opp_decl_dat()` but takes a **file** instead of **data**, reading in the dataset from the HDF5 file using the keyword **name**.

   .. warning::
      The number of data elements specified by the **dim** parameter must match the number of data elements present in the HDF5 file.

