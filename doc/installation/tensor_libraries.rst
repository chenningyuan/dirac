:orphan:

Notes on setting up and using tensor Libraries in Dirac
-------------------------------------------------------


The ExaCorr code (see :ref:`**EXACC`) is built around the use of tensor contract libraries.

Originally based on the tensor libraries `TAL-SH <https://github.com/DmitryLyakh/TAL_SH>`_ and `ExaTENSOR <https://github.com/ORNL-QCI/ExaTENSOR>`_ by Dmitry Lyakh (which support both CPU and GPU architectures), through the Tensor Algebra Processing Primitives (`TAPP <https://arxiv.org/abs/2601.07827>`__) API the code now supports additional tensor libraries: `TBLIS <https://github.com/MatthewsResearchGroup/tblis>`__ (CPU architectures), and the `TAPP reference implementation <https://github.com/TAPPorg/reference-implementation>`__ (CPU architectures, testing only). 

The choice of tensor library used as computational backend is made at compile time, via the ``TENSOR_EXECUTOR`` environment variable.

Configuring for TBLIS (build will fetch TAPP and TBLIS from their repositories)::

    export TENSOR_EXECUTOR=1 # enables TBLIS via the TAPP API

The variable ``ENABLE_TBLIS`` can be toggled (=ON/OFF) through CMake. 

Configuring for the TAPP reference implementation (build will fetch TAPP and the reference implementation from its repository)::

    export TENSOR_EXECUTOR=4 # enables the TAPP reference implementation via the TAPP API

Configuring for TAL-SH (build will fetch the code from the ExaTENSOR repository)::

    export TENSOR_EXECUTOR=2 # enables TAL-SH 

Configuring for ExaTENSOR (build will fetch the code from its repository)::

    export TENSOR_EXECUTOR=3 # enables ExaTENSOR 

For TAL-SH and ExaTENSOR, it is possible to avoid fetching via the network by setting the ``EXATENSOR_GIT_REPO_LOCATION`` variable on CMake to the path to a local clone of the  
ExaTENSOR git repository (the `fork maintained by the DIRAC team <https://github.com/RelMBdev/ExaTENSOR>`__  should be used). This is aimed as a workaround in systems which do not
allow network access, such as some supercomputer centers.

In order to enable GPU support for TAL-SH or ExaTENSOR, the environment variable ``EXA_GPUS`` should be set. 

For NVIDIA GPUs::

    export EXA_GPUS=NVIDIA

Upon starting the build for TAL-SH or ExaTENSOR, the file  ``Exatensor_ENV`` is created, containing the configuration guessed by CMake. It is particularly important, for NVIDIA architectures, to see whether the ``GPU_SM_ARCH`` is correctly set. 

For AMD GPUs::

    export EXA_GPUS=AMD

For AMD GPUs, ``Exatensor_ENV`` will contain ``GPU_CUDA=CUDA`` and ``USE_HIP=YES``. The ``GPU_SM_ARCH`` will be disregarded if present, and the build system should detect the target architecture from environment variables; if not, the GPU type can be passed on via the environment variable ``HIP_TARGET``. 


Running ExaCorr in parallel and on GPU nodes
--------------------------------------------

DIRAC can be compiled and run in parallel with all of the above tensor libraries, provided that the code is compile with 32-bit integers.  However, unless the ExaTENSOR library is use to distribute tensors over multiple nodes, the tensor libraries will only use a single node.

In the case ExaCorr has been compiled with GPU support for TAL-SH, the environment variable ``TALSH_GPUS`` can be used to control how many of the GPUs accessible to TAL-SH runs will be used (if ``TALSH_GPUS`` is not defined, all GPUs available to the process will be used).

ExaTENSOR runs will require the definition of another set of environment variables, for example as in::

    export QF_NUM_PROCS=4             # total number of MPI ranks in the calculation. Must be at least 4 for ExaTENSOR to work properly. 
    export QF_PROCS_PER_NODE=4        # number of MPI ranks executing per node
    export QF_CORES_PER_PROCESS=2     # Number of CPU cores attached to each MPI rank
    export QF_NVMEM_PER_PROCESS=0     #
    export QF_HOST_BUFFER_SIZE=1000   # Memory available on the node, divided by the number of MPI ranks
    export QF_MEM_PER_PROCESS=750     # about 75% of QF_HOST_BUFFER_SIZE

    export QF_GPUS_PER_PROCESS=n # Number of NVIDIA or AMD GPUs to be used per MPI rank
    export QF_MICS_PER_PROCESS=0 # this should always be set to zero
    export QF_AMDS_PER_PROCESS=0 # this should always be set to zero, even in the case of AMD GPU systems

    export OMP_NUM_THREADS=$QF_CORES_PER_PROCESS
    export OMP_DYNAMIC=false
    export OMP_MAX_ACTIVE_LEVELS=3
    export OMP_THREAD_LIMIT=256
    export OMP_WAIT_POLICY=PASSIVE
    export OMP_PROC_BIND="spread"
    export OMP_PLACES="{0:4},{4:4},{8:4},{12:4}"

 (please consult the examples for runtime configurations in the source tree, and the documentation of ExaTENSOR for further details)



