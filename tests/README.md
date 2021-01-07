The test case is a 3D mesh generated via distributed computing.

### How to test Microsoft MPI (https://docs.microsoft.com/en-us/message-passing-interface/microsoft-mpi)

#### go to the Solution Directory(where sln file is) > WinDolfin\3rdparties\MSMPI\Bin
#### Run " .\mpiexec.exe -n [Number of Processes] Solution Directory\WinDolfin\x64\Debug\tests.exe"


##### MPI_Init_thread level = MPI_THREAD_MULTIPLE
##### main is_receiver 2 Build LocalMeshData from local Mesh
##### MPI_Init_thread level = MPI_THREAD_MULTIPLE
##### main is_receiver 1 Build LocalMeshData from local Mesh
##### sub is_receiver 2
##### receive_mesh_data
##### sub is_receiver 1
##### receive_mesh_data
##### MPI_Init_thread level = MPI_THREAD_MULTIPLE
##### main is_broadcaster 0 Build LocalMeshData from local Mesh
##### sub is_broadcaster 0
##### extract_mesh_data geometry.dim2
##### broadcast_mesh_data
##### Received %d vertex coordinates 5547
##### Received %d vertex coordinates 5547
##### Received %d vertex coordinates 5547
##### Received %d cell coordinates 5461
##### Received %d cell coordinates 5462
##### Attempting to use an MPI routine after finalizing MPI

### How to install MSMPI
#### Download (https://www.microsoft.com/en-us/download/details.aspx?id=100593)
#### Install both msmpisetup.exe and msmpisdk.msi
#### Copy the Installed Directory\Microsoft MPI\Bin to Solution Directory\WinDolfin\3rdparties\MSMPI