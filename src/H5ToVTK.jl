module H5ToVTK

using HDF5
using WriteVTK

export h5tovtk

function convert_cell_type(s::String)
    if s == "tetrahedron"
        return VTKCellTypes.VTK_TETRA
    elseif s == "triangle"
        return VTKCellTypes.VTK_TRIANGLE
    else
        @error "Unknown cell type: $s"
    end
end

function h5tovtk(filename::String, groups::Vector{String}; coords_tag::String="coordinates", topo_tag::String="topology")
    fid = h5open(filename * ".h5", "r")
    vtm = vtk_multiblock(filename)
    for grp in groups
        block = multiblock_add_block(vtm, filename * "_" * grp)
        h5_grp = fid[grp]
        coords = read(h5_grp, coords_tag) 
        topo_obj = h5_grp[topo_tag]
        single_cell_type = convert_cell_type(read(topo_obj["celltype"]))
        cells = map(eachcol(read(topo_obj))) do v
            MeshCell(single_cell_type, v .+ 1)
        end
        if haskey(h5_grp, "values")
            vals = read(h5_grp, "values")
            map(unique(vals)) do dif_v 
                ids = findall(==(dif_v), vals)
                vtk_grid(block, coords, cells[ids])
            end
        else 
            vtk_grid(block, coords, cells)
        end
    end
    close(vtm)
    return nothing
end

end
