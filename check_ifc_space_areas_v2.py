#!/usr/bin/env python3
"""Compare IfcSpace area attributes with a geometry-derived horizontal floor-face area.


Usage:
python check_ifc_space_areas_v2.py FILENAME.ifc --threshold-mm2 500


Requires IfcOpenShell, NumPy, pandas, openpyxl.

Method:
- Uses IFC geometry in world coordinates.
- Triangulates the mesh and identifies nearly-horizontal triangles (normal Z close to 1).
- Clusters those triangles by elevation and uses the lowest horizontal-face elevation
  cluster as the floor candidate (not all room surfaces / not a closed-mesh /2 trick).
- Computes the candidate floor area as the sum of projected triangle areas.
- Reports only differences above threshold, plus missing attributes/geometry issues.

IMPORTANT: The method is still a geometric heuristic. IFC space meshes are not
uniformly exported with a dedicated floor face. Inspect the "Geometry diagnostics"
worksheet for floor candidates and confidence notes before treating flags as errors.
"""
from __future__ import annotations

import argparse
from pathlib import Path

from collections import defaultdict
import numpy as np
import pandas as pd

import ifcopenshell
import ifcopenshell.geom


AREA_NAMES = ("NetFloorArea", "GrossFloorArea", "Area", "RoomArea", "SpaceArea")
HORIZONTAL_NORMAL_Z = 0.985
DEFAULT_ELEVATION_TOL_M = 0.015





def read_area_properties(space):
    """Return preferred area property/quantity and all numeric area-like values."""
    values = {}
    for rel in getattr(space, "IsDefinedBy", []) or []:
        definition = getattr(rel, "RelatingPropertyDefinition", None)
        if not definition:
            continue
        if definition.is_a("IfcPropertySet"):
            for prop in getattr(definition, "HasProperties", []) or []:
                if not prop.is_a("IfcPropertySingleValue"):
                    continue
                wrapped = getattr(getattr(prop, "NominalValue", None), "wrappedValue", None)
                if isinstance(wrapped, (int, float)) and not isinstance(wrapped, bool):
                    name = str(getattr(prop, "Name", ""))
                    if "area" in name.lower() or "floor" in name.lower() or "площад" in name.lower():
                        values[name] = (float(wrapped), str(getattr(definition, "Name", "")))
        elif definition.is_a("IfcElementQuantity"):
            for quantity in getattr(definition, "Quantities", []) or []:
                if quantity.is_a("IfcQuantityArea"):
                    values[str(quantity.Name)] = (float(quantity.AreaValue), str(definition.Name))

    for preferred in AREA_NAMES:
        if preferred in values:
            return preferred, values[preferred][0], values[preferred][1], values
    for name, (value, pset) in values.items():
        if "area" in name.lower() or "floor" in name.lower() or "площад" in name.lower():
            return name, value, pset, values
    return None, None, None, values



def horizontal_floor_candidate(space, settings, elevation_tol_m=DEFAULT_ELEVATION_TOL_M):
    """Return area of lowest cluster of horizontal triangles, plus diagnostics.

    The mesh is expected to be in metres (IfcOpenShell geometry convention).
    Triangle area is projected on XY; sloped surfaces are excluded by normal filter.
    """
    shape = ifcopenshell.geom.create_shape(settings, space)
    geom = shape.geometry
    vertices = np.asarray(geom.verts, dtype=np.float64).reshape((-1, 3))
    faces = np.asarray(geom.faces, dtype=np.int64).reshape((-1, 3))
    if len(faces) == 0 or len(vertices) == 0:
        return None, {"Status": "No triangles", "Horizontal triangles": 0}


    p0, p1, p2 = vertices[faces[:, 0]], vertices[faces[:, 1]], vertices[faces[:, 2]]
    cross = np.cross(p1 - p0, p2 - p0)
    twice_area = np.linalg.norm(cross, axis=1)
    valid = twice_area > 1e-12
    nz_ratio = np.zeros_like(twice_area)
    nz_ratio[valid] = np.abs(cross[valid, 2]) / twice_area[valid]
    horiz = valid & (nz_ratio >= HORIZONTAL_NORMAL_Z)
    if not np.any(horiz):
        return None, {"Status": "No near-horizontal faces", "Horizontal triangles": 0,
                     "Mesh triangles": len(faces)}

    # Each triangle's centroid elevation. Faces at nearly the same z are grouped.
    cent_z = ((p0[:, 2] + p1[:, 2] + p2[:, 2]) / 3.0)[horiz]
    proj_area = (np.abs(cross[:, 2]) / 2.0)[horiz]
    order = np.argsort(cent_z)
    cent_z, proj_area = cent_z[order], proj_area[order]

    # Cluster sorted elevations using a fixed tolerance, updating cluster mean.
    clusters = []
    current_z, current_areas = [], []
    for z, area in zip(cent_z, proj_area):
        if not current_z or abs(z - np.mean(current_z)) <= elevation_tol_m:
            current_z.append(float(z)); current_areas.append(float(area))
        else:
            clusters.append((float(np.mean(current_z)), float(sum(current_areas)), len(current_z)))
            current_z, current_areas = [float(z)], [float(area)]
    if current_z:
        clusters.append((float(np.mean(current_z)), float(sum(current_areas)), len(current_z)))

    # Lowest horizontal cluster is the floor candidate; unlike prior code, no division by 2.
    floor_z, floor_area, floor_triangle_count = clusters[0]
    diagnostic = {
        "Status": "OK - lowest horizontal cluster used",
        "Floor candidate elevation (m)": floor_z,
        "Floor candidate area (m²)": floor_area,
        "Horizontal triangles": int(horiz.sum()),
        "Floor triangles": floor_triangle_count,
        "Horizontal elevation clusters": len(clusters),
        "Mesh triangles": len(faces),
        "All horizontal clusters (z m : area m²)": " | ".join(f"{z:.4f}: {a:.6f} ({n} tri)" for z, a, n in clusters),
    }
    return floor_area, diagnostic



def main():
    parser = argparse.ArgumentParser(description="Check IfcSpace area properties against horizontal floor-face geometry.")
    parser.add_argument("ifc", help="Input IFC file")
    parser.add_argument("--threshold-mm2", type=float, default=500.0,
                       help="Report only area differences greater than this value (default 500 mm²)")
    parser.add_argument("--elevation-tolerance-mm", type=float, default=15.0,
                       help="Cluster horizontal triangles by elevation (default 15 mm)")
    parser.add_argument("--output", default="ifc_space_area_mismatches_v2.xlsx", help="Output Excel path")
    args = parser.parse_args()



    threshold_m2 = args.threshold_mm2 / 1_000_000.0
    elevation_tol_m = args.elevation_tolerance_mm / 1000.0
    model = ifcopenshell.open(args.ifc)
    spaces = model.by_type("IfcSpace")
    settings = ifcopenshell.geom.settings()
    settings.set(settings.USE_WORLD_COORDS, True)


    mismatches, diagnostics, errors = [], [], []
    for space in spaces:
        gid = getattr(space, "GlobalId", "")
        name = getattr(space, "Name", "")
        long_name = getattr(space, "LongName", "")
        attr_name, attr_area, pset, all_values = read_area_properties(space)
        try:
            geom_area, diag = horizontal_floor_candidate(space, settings, elevation_tol_m)
        except Exception as exc:
            geom_area = None
            diag = {"Status": "Geometry exception", "Error": str(exc)}
            errors.append({"ExpressID": space.id(), "GlobalId": gid, "Name": name, "Error": str(exc)})

        diag.update({"ExpressID": space.id(), "GlobalId": gid, "Name": name, "Long Name": long_name})
        diagnostics.append(diag)

        if attr_area is not None and geom_area is not None:
            difference = abs(attr_area - geom_area)
            if difference > threshold_m2:
                mismatches.append({
                    "ExpressID": space.id(), "GlobalId": gid, "Name": name, "Long Name": long_name,
                    "Attribute": attr_name, "Property Set": pset,
                    "Attribute Area (m²)": attr_area, "Floor Candidate Area (m²)": geom_area,
                    "Difference (m²)": difference, "Difference (mm²)": difference * 1_000_000,
                    "Other area properties": "; ".join(f"{k}={v[0]} ({v[1]})" for k, v in all_values.items()),
                    "Geometry status": diag.get("Status", ""),
                })
        elif attr_area is None:
            mismatches.append({
                "ExpressID": space.id(), "GlobalId": gid, "Name": name, "Long Name": long_name,
                "Status": "Area property not found", "Floor Candidate Area (m²)": geom_area,
            })
        elif geom_area is None:
            mismatches.append({
                "ExpressID": space.id(), "GlobalId": gid, "Name": name, "Long Name": long_name,
                "Status": "Floor candidate could not be calculated", "Attribute": attr_name,
                "Attribute Area (m²)": attr_area, "Geometry status": diag.get("Status", ""),
            })

    summary = pd.DataFrame([{
        "IFC file": str(Path(args.ifc).resolve()),
        "IfcSpace count": len(spaces),
        "Threshold (mm²)": args.threshold_mm2,
        "Threshold (m²)": threshold_m2,
        "Elevation cluster tolerance (mm)": args.elevation_tolerance_mm,
        "Mismatch / missing rows": len(mismatches),
        "Geometry exceptions": len(errors),
        "Method": "Sum projected area of near-horizontal triangles at the lowest elevation cluster; no mesh-wide division by 2.",
        "Limitation": "Lowest horizontal mesh faces may not be the room floor in every IFC export. Review Geometry Diagnostics before assuming a flag is an error.",
    }])

    with pd.ExcelWriter(args.output, engine="openpyxl") as writer:
        summary.to_excel(writer, index=False, sheet_name="Summary")
        pd.DataFrame(mismatches).to_excel(writer, index=False, sheet_name="Mismatches")
        pd.DataFrame(diagnostics).to_excel(writer, index=False, sheet_name="Geometry Diagnostics")
        pd.DataFrame(errors).to_excel(writer, index=False, sheet_name="Geometry Errors")
        # Format columns, freeze top row, and auto-filter
        for ws in writer.book.worksheets:
            ws.freeze_panes = "A2"
            ws.auto_filter.ref = ws.dimensions
            for column_cells in ws.columns:
                col_letter = column_cells[0].column_letter
                max_len = min(60, max((len(str(c.value)) for c in column_cells if c.value is not None), default=10))
                ws.column_dimensions[col_letter].width = max(12, max_len + 2)

    print(f"IfcSpace: {len(spaces)}")
    print(f"Rows in mismatches report: {len(mismatches)}")
    print(f"Geometry exceptions: {len(errors)}")
    print(f"Report written: {Path(args.output).resolve()}")

if __name__ == "__main__":
    main()
