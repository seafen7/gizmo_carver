"""
   writer_gizmo_carver.py

   Purpose:
        Writer class that compresses the relevant layers needed to create  
        RADMC amr, number density, line, and dust files. Should not need to edit
        this file.

   Author:
        Sean Feng, feng.sean01@utexas.edu
        Spring 2022
        
        Modified from: carver_CarveOut.py, written by:
        Aaron T. Lee, aaron.t.lee@utexas.edu
        Spring 2018

   Written/Tested with Python 3.9, yt 4.0.2
"""

import numpy as np
import yt
from yt.utilities.lib.write_array import write_3D_array, write_3D_vector_array # same as carver_CarveOut
from yt.extensions.astro_analysis.radmc3d_export.RadMC3DInterface import RadMC3DLayer, RadMC3DSource

class RadMC3DWriter_Gizmo:
    
    # From yt GitHub RadMC3DInterface.py RadMC3dWriter class
    def __init__(self, ds, a_boxLeft=None, a_boxRight=None, a_boxDim=0):
        self.max_level = 0
        self.cell_count = 0
        self.layers = []
        self.domain_dimensions = [int(np.abs(a_boxDim)), int(np.abs(a_boxDim)), int(np.abs(a_boxDim))]# ds.domain_dimensions
        self.domain_left_edge = ds.domain_left_edge
        self.domain_right_edge = ds.domain_right_edge
        self.ds = ds
        self.force_unigrid = 1
        
        if a_boxLeft is None:
            self.boxLeft = ds.domain_left_edge
        else:
            self.boxLeft = a_boxLeft
        
        if a_boxRight is None:
            self.boxRight = ds.domain_right_edge
        else:
            self.boxRight = a_boxRight
            
        if a_boxDim == 0:
            self.boxDim = 256j
        else:
            self.boxDim = a_boxDim
        
        # Code to add more layers for higher resolution
        # if(self.force_unigrid):
        #     #self.max_level = 0
        #     newDim      = [(2**a_max_level)*x for x in a_boxDim]  # Need to specify Ncells and box_dim
        #     self.domain_dimensions = newDim
        #     base_layer = RadMC3DLayer(0, None, 0,
        #                               self.domain_left_edge,
        #                               self.domain_right_edge,
        #                               newDim,
        #                               self.allow_periodic, self.domainPatches)
        #     self.cell_count += np.product(newDim) #np.product(a_ds.domain_dimensions)

        base_layer = RadMC3DLayer(
            0,
            None,
            0,
            self.domain_left_edge,
            self.domain_right_edge,
            self.domain_dimensions,
        )

        self.layers.append(base_layer)
        self.cell_count += np.product(self.domain_dimensions)


        # sorted_grids = sorted(ds.index.grids, key=lambda x: x.Level)
        # for grid in sorted_grids:
        #     if grid.Level <= self.max_level:
        #         self._add_grid_to_layers(grid)
        
        
    # From carver_CarveOut.py
    def write_amr_grid(self, filename):
        '''
        This routine writes the "amr_grid.inp" file that describes the mesh
        radmc3d will use.
    
        '''
        dims = self.domain_dimensions # if force_unigrid, this should have been overwritten to the right value
    
        LE = self.boxLeft # carved out region
        RE = self.boxRight
    
        #CellCount = np.product(self.domain_dimensions) # check
    
        # Taken from YT, fairly certain shouldn't be necessary, since O2 code_length is CGS
        # RadMC-3D wants the cell wall positions in cgs. Convert here:
        LE_cgs = LE.in_units('cm').d  # don't write the units, though
        RE_cgs = RE.in_units('cm').d  # also prevents passing by pointer (self.domain... etc. needs to be used again)
    
        # Shift so centered at 0,0,0
        Center_cgs = 0.5*(LE_cgs + RE_cgs)
        LE_cgs = LE_cgs - Center_cgs
        RE_cgs = RE_cgs - Center_cgs
    
        # calculate cell wall positions (may potentially exceed full domain if periodic)
        xs = [str(x) for x in np.linspace(LE_cgs[0], RE_cgs[0], dims[0]+1)]
        ys = [str(y) for y in np.linspace(LE_cgs[1], RE_cgs[1], dims[1]+1)]
        zs = [str(z) for z in np.linspace(LE_cgs[2], RE_cgs[2], dims[2]+1)]
    
        # writer file header
        grid_file = open(filename, 'w')
        grid_file.write('1 \n')  # iformat is always 1
        if(self.max_level == 0 or self.force_unigrid==1):
            grid_file.write('0 \n')
        else:
            grid_file.write('10 \n')  # only layer-style files are supported
        grid_file.write('0 \n') # '1 \n') # only cartesian coordinates are supported
        grid_file.write('0 \n')
        grid_file.write('{}    {}    {} \n'.format(1, 1, 1))  # assume 3D
        grid_file.write('{}    {}    {} \n'.format(dims[0], dims[1], dims[2]))
        if(self.max_level != 0 and self.force_unigrid==0):
            s = str(self.max_level) + '    ' + str(len(self.layers)-1) + '\n'
            grid_file.write(s)
    
        # write base grid cell wall positions
        for x in xs:
            grid_file.write(x + '    ')
        grid_file.write('\n')
    
        for y in ys:
            grid_file.write(y + '    ')
        grid_file.write('\n')
    
        for z in zs:
            grid_file.write(z + '    ')
        grid_file.write('\n')
    
        # write information about fine layers
        for layer in self.layers[1:]: # [1:] skips entry 0, the base layer
            p = layer.parent # parent id
            dds = (layer.RightEdge - layer.LeftEdge) / (layer.ActiveDimensions) # cell size
            if p == 0: # parent is the base layer,
                # have to incorporate periodicity, if used
                # LE = base layer, needs to be adjusted
                # layer.Left is the level=1 layer, no adjustment needed
                # Beginning index
                ind = (layer.LeftEdge - LE) / (2.0*dds) + 1 # LE wasn't shifted so 0,0,0 is center. OK!
                if(layer.is_periodic):
                    for shifted_grid in self.domainPatches:
                        ledge_shift, redge_shift = shifted_grid
                        LE = np.maximum(yt.YTArray(ledge_shift,'cm'),  layer.LeftEdge)
                        RE = np.minimum(yt.YTArray(redge_shift,'cm'), layer.RightEdge)
                        if np.any(RE > LE):
                            print("HOW COULD THIS HAPPEN?!?!")
                            LE = yt.YTArray(ledge_shift,'cm')
                            ind = (layer.LeftEdge - LE) / (2.0*dds) + 1
                            #print(ind)
                            break # ASSUMPTION: it can overlap with only one of the broken up boxes
            else: # parent is an AMR layer
                parent_LE = np.zeros(3)
                for potential_parent in self.layers:
                    if potential_parent.id == p:
                        parent_LE = potential_parent.LeftEdge
                        break # might save a few tick tocks
                ind = (layer.LeftEdge - parent_LE) / (2.0*dds) + 1 # #index in parent grid cells; periodic or not, they both should be correct; difference is good here
            #print(str(np.array(ind)+0.5) + " " + str([int(x) for x in np.array(ind)+0.5]))
            ix = int(ind[0]+0.5) # beginning index in terms of parent grid cells (b/c of the 2*dds)
            iy = int(ind[1]+0.5)
            iz = int(ind[2]+0.5)
            #ix = math.ceil(ind[0]) # beginning index in terms of parent grid cells (b/c of the 2*dds)
            #iy = math.ceil(ind[1])
            #iz = math.ceil(ind[2])
            nx, ny, nz = layer.ActiveDimensions / 2 # number of cells (/2 to measure in dimensions of parent cell size)
            #print(nx,ny,nz)
            s = '{}    {}    {}    {}    {}    {}    {} \n'
            #s = s.format(p, ix, iy, iz, nx, ny, nz)
            #s = s.format(p, ix, iy, iz, int(np.rint(nx.d)), int(np.rint(ny.d)), int(np.rint(nz.d))) # RAD complains if anything is not int()
            s = s.format(p, ix, iy, iz, int(nx), int(ny), int(nz) ) # RAD complains if anything is not int()
            #s = s.format(p, ix, iy, iz, int(round(nx)), int(round(ny)), int(round(nz)) ) # RAD complains if anything is not int()
            #CellCount = CellCount + 2*(round(nx)*round(ny)*round(nz))
            grid_file.write(s)
    
        grid_file.close()
        
        
    # Custom covering_grid function
    def _covering_grid(self, ds):
        le = self.boxLeft
        re = self.boxRight
        res = self.boxDim
        
        grid = ds.r[le[0]:re[0]:res, le[1]:re[1]:res, le[2]:re[2]:res]
        
        return grid
    
    
    # From yt GitHub RadMC3DInterface.py 
    def _write_layer_data_to_file(self, fhandle, field, level, LE, dim):
        cg = self._covering_grid(self.ds) # Replaced with custom covering_grid function
        
        if isinstance(field, list):
            data_x = cg[field[0]]
            data_y = cg[field[1]]
            data_z = cg[field[2]]
            write_3D_vector_array(data_x, data_y, data_z, fhandle)
        else:
            data = cg[field]
            write_3D_array(data, fhandle)
                
    
    # From yt GitHub RadMC3DInterface.py 
    def write_line_file(self, field, filename):
            """
            This method writes out fields in the format radmc3d needs to compute
            line emission.
            Parameters
            ----------
            field : string or list of 3 strings
                If a string, the name of the field to be written out. If a list,
                three fields that will be written to the file as a vector quantity.
            filename : string
                The name of the file to write the data to. The filenames radmc3d
                expects for its various modes of operation are described in the
                radmc3d manual.
            """
            fhandle = open(filename, "w")
    
            # write header
            fhandle.write("1 \n")
            fhandle.write(str(self.cell_count) + " \n")
    
            # now write fine layers:
            for layer in self.layers:
                lev = layer.level
                if lev == 0:
                    LE = self.boxLeft # domain_left_edge
                    N = self.domain_dimensions
                else:
                    LE = layer.LeftEdge
                    N = layer.ActiveDimensions
    
                self._write_layer_data_to_file(fhandle, field, lev, LE, N)
    
            fhandle.close()
        
        
    # From yt GitHub RadMC3DInterface.py 
    def write_dust_file(self, field, filename):
        """
        This method writes out fields in the format radmc3d needs to compute
        thermal dust emission. In particular, if you have a field called
        "DustDensity", you can write out a dust_density.inp file.
        Parameters
        ----------
        field : string
            The name of the field to be written out
        filename : string
            The name of the file to write the data to. The filenames radmc3d
            expects for its various modes of operations are described in the
            radmc3d manual.
        """
        fhandle = open(filename, "w")

        # write header
        fhandle.write("1 \n")
        fhandle.write(str(self.cell_count) + " \n")
        fhandle.write("1 \n")

        # now write fine layers:
        for layer in self.layers:
            lev = layer.level
            if lev == 0:
                LE = self.boxLeft # domain_left_edge
                N = self.domain_dimensions
            else:
                LE = layer.LeftEdge
                N = layer.ActiveDimensions

            self._write_layer_data_to_file(fhandle, field, lev, LE, N)

        fhandle.close()
        
        
        
        