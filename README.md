# thorax_iphone_to_FEMmesh
Construct a subject specific FEM mesh with encoded electrodes from a segmented/processed 3D iPhone scan of a subject's thorax. The main function, construct_FEMmeshes.m, 
  *  Loads a processed 3D scan, processed via https://github.com/ethankmurphy/EIT_3D_thorax_scans 
  *  Then it construct an outer surface triangulation with encoded electrodes.
      * The meshing here is done using distmesh, http://persson.berkeley.edu/distmesh/#:~:text=DistMesh%20is%20a%20simple%20MATLAB,Department%20of%20Mathematics%20at%20MIT.
      * One needs just to download distmesh and add the path (lines 28-31 of construct_FEMmeshes.m)
  * After the surface mesh is made, then the information is sent to gmsh, where a 3D tetrahedral mesh is generated.
      * gmsh is a free software, https://gmsh.info/
      * One needs to download/install the software and then supply the code the path to the executable (see lines 22-26 in construct_FEMmeshes.m)
      
The final mesh should look something like this (see below).  
    
![example_processed_FEMmesh](https://user-images.githubusercontent.com/87721360/219524814-b670648a-cefa-4d75-bc5c-7383d0aff8e3.png)
