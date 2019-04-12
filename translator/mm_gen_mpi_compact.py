# Open source copyright declaration based on BSD open source template:
# http://www.opensource.org/licenses/bsd-license.php
#
# This file is part of the OPS distribution.
#
# Copyright (c) 2013, Mike Giles and others. Please see the AUTHORS file in
# the main source directory for a full list of copyright holders.
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
# Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
# The name of Mike Giles may not be used to endorse or promote products
# derived from this software without specific prior written permission.
# THIS SOFTWARE IS PROVIDED BY Mike Giles ''AS IS'' AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL Mike Giles BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

"""
MM MPI_seq code generator (for C/C++ applications)

This routine is called by ops.py which parses the input files

It produces a file xxx_seq_kernel.cpp for each kernel,
plus a master kernel file

"""

import re
import datetime
import os
import glob

import util
import config

para_parse = util.para_parse
comment_remover = util.comment_remover
remove_trailing_w_space = util.remove_trailing_w_space
parse_signature = util.parse_signature
check_accs = util.check_accs
mult = util.mult

comm = util.comm
code = util.code
FOR = util.FOR
FOR2 = util.FOR2
WHILE = util.WHILE
ENDWHILE = util.ENDWHILE
ENDFOR = util.ENDFOR
IF = util.IF
ELSEIF = util.ELSEIF
ELSE = util.ELSE
ENDIF = util.ENDIF


def mm_gen_mpi_compact(master, date, kernels):

  GBL   = 2;

  READ = 1;  WRITE = 2;  RW  = 3;
  INC  = 4;  MAX   = 5;  MIN = 6;

  accsstring = ['READ','WRITE','RW','INC','MAX','MIN' ]

##########################################################################
#  create new kernel file
##########################################################################

  for nk in range (0,len(kernels)):
    arg_type  = kernels[nk]['arg_type']
    arg_type2  = kernels[nk]['arg_type2']
    arg_qtype  = kernels[nk]['arg_qtypes']
    name  = kernels[nk]['name']
    nargs = kernels[nk]['nargs']
    dim   = kernels[nk]['dim']
    dims  = kernels[nk]['dims']
    stens = kernels[nk]['stens']
    var   = kernels[nk]['var']
    accs  = kernels[nk]['accs']
    spaces  = kernels[nk]['spaces']
    space  = kernels[nk]['space']
    lam   = kernels[nk]['lam']

    config.file_text = ''
    config.depth = 0
    n_per_line = 4
    space_counters=['i','j','k','l','m','n']

##########################################################################
#  start with seq kernel function
##########################################################################
    code('#include <iostream>\n#include "MultidimArray.hpp"\n#include "compressed_cell_centric/Arguments.hpp"\n#include "compressed_cell_centric/Computation.hpp"\n\nusing namespace MM;\nusing namespace MM::compressed_cell_centric;')
    code('')
    code('#define NeighProxy NeighProxyDirect')
    code('')
    comm('user function')
    pos = lam.find('(')
    code('inline void '+name+lam[pos:])

    j = lam.find('{',pos)
    k = para_parse(lam, j, '{', '}')
    arg_list = parse_signature(lam[pos:j])
    print arg_list

    code('')

##########################################################################
#  now host stub
##########################################################################

    comm(' host stub function')
    code('void mm_par_loop_'+name+'(std::string name, Computation<'+str(dim)+'> &computation,')
    text = ''
    for n in range (0, nargs):
      text = text +' '+arg_qtype[n]+' arg'+str(n)
      if nargs <> 1 and n != nargs-1:
        text = text +','
      else:
        text = text +') {'
      if n%n_per_line == 3 and n <> nargs-1:
         text = text +'\n'
    code(text);
    config.depth = 2

    code('const std::array<std::size_t, '+str(dim)+'> &begin = computation.index_generator.get_begin();')
    code('const std::array<std::size_t, '+str(dim)+'> &end = computation.index_generator.get_end();')
    code('')
    code('CompressedDataStructure &structure = computation.data.structure;')
    code('std::array<std::size_t,'+str(dims[n])+'> shape = computation.data.get_size();')

    for n in range(0,nargs):
        if arg_type[n] == 'arg_dat':
            code('double * __restrict data'+str(n)+' = arg'+str(n)+'.get_raw();')
            if 'CellMatData' == spaces[n]:
              code('double *__restrict data'+str(n)+'_list = arg'+str(n)+'.get_raw_list();')

    if 'Cell' in space:
      line = ''
      code('#pragma omp parallel for collapse('+str(dim-1)+')')
      for d2 in range(0,dim):
        d = dim-d2-1
        if d2 == dim-1:
            code('#pragma omp simd')
        FOR('std::size_t '+space_counters[d]+' = begin['+str(d)+']; '+space_counters[d]+' < end['+str(d)+']; '+space_counters[d]+'++')
        line = space_counters[d]+', '+line
    else:
      line = ''
      for d2 in range(0,dim):
        line = line + '0, '
    coords = 'Coords<'+str(dim)+'>('+line[:-2]+')'

    if 'CellMat' in space:
      IF('structure.structure[i+j*shape[0]].nmats == 1') #single MAT
      code('std::size_t mat_index = structure.structure[i+j*shape[0]].imat;')
      #
      # Function call
      #
      for arg in range(0,nargs):
          if arg_type[arg] == 'arg_dat':
              if arg_type2[arg] == 'NEIGH':
                code(arg_list[arg]+'(computation.data, '+coords+', CellMatIndex(i+j*shape[0], mat_index), ValueIndex(ValueIndex::Type::SINGLE_MAT, i+j*shape[0]), arg'+str(arg)+'.dataset);')
              elif 'CellData' == spaces[arg]:
                code(arg_list[arg]+'(data'+str(arg)+'[i+j*shape[0]]);')
              elif 'MatData' == spaces[arg]:
                code(arg_list[arg]+'(data'+str(arg)+'[mat_index]);')
              else:
                code(arg_list[arg]+'(data'+str(arg)+'[i+j*shape[0]]);')
          else:
            code(arg_list[arg]+'(arg'+str(arg)+'.get('+coords+', computation.data, CellMatIndex(i+j*shape[0], mat_index), ValueIndex(ValueIndex::Type::SINGLE_MAT, i+j*shape[0])));')
      code(lam[j+1:k])
      ELSE()

      # Using a linked list in the compressed cell centric mixed storage if MM_LINKED is defined.
      code('#ifdef MM_LINKED')
      FOR('std::size_t structure_index = structure.structure[i+j*shape[0]].imat;\n\
            structure_index!=-1ul; structure_index = structure.mixed_storage[structure_index].nextfrac')
      code('#else')
      code('const Cell& cell = structure.structure[i+j*shape[0]];')
      FOR('std::size_t structure_index = cell.imat;\n\
             structure_index < cell.imat + cell.nmats; structure_index++')
      code('#endif')
      code('std::size_t mat_index = structure.mixed_storage[structure_index].material;')

      #
      # Function call
      #
      for arg in range(0,nargs):
          if arg_type[arg] == 'arg_dat':
              if arg_type2[arg] == 'NEIGH':
                code(arg_list[arg]+'(computation.data, '+coords+', CellMatIndex(i+j*shape[0], mat_index), ValueIndex(ValueIndex::Type::MULTIMAT, structure_index), arg'+str(arg)+'.dataset);')
              elif 'CellData' == spaces[arg]:
                code(arg_list[arg]+'(data'+str(arg)+'[i+j*shape[0]]);')
              elif 'MatData' == spaces[arg]:
                code(arg_list[arg]+'(data'+str(arg)+'[mat_index]);')
              else:
                code(arg_list[arg]+'(data'+str(arg)+'_list[structure_index]);')
          else:
            code(arg_list[arg]+'(arg'+str(arg)+'.get('+coords+', computation.data, CellMatIndex(i+j*shape[0], mat_index), ValueIndex(ValueIndex::Type::MULTIMAT, structure_index)));')
      code(lam[j+1:k])
      ENDFOR()
      ENDIF()
    elif 'Mat' in space:
      FOR('std::size_t mat_index = 0; mat_index < computation.data.get_mat_number(); mat_index++')
    else:
      code('std::size_t mat_index = 0;')



    if 'Cell' in space:
      for d2 in range(0,dim):
        ENDFOR()

    config.depth = 0
    code('}')

##########################################################################
#  output individual kernel file
##########################################################################
    if not os.path.exists('./MPI'):
      os.makedirs('./MPI')
    fid = open('./MPI/'+name+'_seq_kernel.cpp','w')
    date = datetime.datetime.now()
    fid.write('//\n// auto-generated by mm.py\n//\n')
    fid.write(config.file_text)
    fid.close()

# end of main kernel call loop
