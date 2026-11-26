# DGP_RoomChekerIFC

Проверка соответсвия Атрибута площади и геометрии IFC помещений 

Compare IfcSpace area attributes with a geometry-derived horizontal floor-face area.


# **Usage**
python check_ifc_space_areas_v2.py FILENAME.ifc --threshold-mm2 500

Requires IfcOpenShell, NumPy, pandas, openpyxl.

# **Method**
- Uses IFC geometry in world coordinates.
- Triangulates the mesh and identifies nearly-horizontal triangles (normal Z close to 1).
- Clusters those triangles by elevation and uses the lowest horizontal-face elevation cluster as the floor candidate (not all room surfaces / not a closed-mesh /2 trick).
- Computes the candidate floor area as the sum of projected triangle areas.
- Reports only differences above threshold, plus missing attributes/geometry issues.

**IMPORTANT**: The method is still a geometric heuristic. IFC space meshes are not
uniformly exported with a dedicated floor face. Inspect the "Geometry diagnostics"
worksheet for floor candidates and confidence notes before treating flags as errors.