# swash 3.0.0

## Breaking changes (Non-backwards compatible)
- Complete change to the nbmatrix() function: creation of an instance of the new nbmatrix class 
- Replacement of the old nbstat() function with the nbstat() method of the nbmatrix class
- Rearrangement of the code in modules for better readability and easier maintaining

## New features
- Spatial statistics based on neighborhood matrix: Global Getis-Ord, Global Moran's I, 
Local Getis-Ord Gi*, with all of them being methods of class nbmatrix
- Plotting a map from an nbmatrix object
- Importing geodata (sf) in infpan objects
- Plotting maps of attributes in an infpan object with method plot_map()

## Bugfixes
- Corrections in RD documentations