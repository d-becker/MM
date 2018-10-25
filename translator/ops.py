#!/usr/bin/env python

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
OPS source code transformation tool (for the C/C++ API)

This tool parses the user's original source code to produce
target-specific code to execute the user's kernel functions.

This prototype is written in Python

usage: ./ops.py file1, file2 ,...

This takes as input

file1.cpp, file2.cpp, ...

and produces as output modified versions

file1_ops.cpp, file2_ops.cpp, ...

then calls a number of target-specific code generators
to produce individual kernel files of the form

xxx_seq_kernel.cpp -- for single threaded x86 execution (also used for MPI)
xxx_omp_kernel.cpp -- for OpenMP x86 execution
xxx_kernel.cu -- for CUDA execution

"""

import sys
import re
import datetime

import util

from mm_gen_mpi import mm_gen_mpi

arithmetic_regex_pattern = r'^[ \(\)\+\-\*\\\.\%0-9]+$'

comment_remover = util.comment_remover
remove_trailing_w_space = util.remove_trailing_w_space

def arg_parse(text, j):
    """Parsing arguments in op_par_loop to find the correct closing brace"""

    depth = 0
    loc2 = j
    arglist = []
    prev_start = j
    while 1:
        if text[loc2] == '(':
            if depth == 0:
                prev_start = loc2+1
            depth = depth + 1

        elif text[loc2] == ')':
            depth = depth - 1
            if depth == 0:
                arglist.append(text[prev_start:loc2].strip())
                return arglist

        elif text[loc2] == ',':
            if depth == 1:
                arglist.append(text[prev_start:loc2].strip())
                prev_start = loc2+1
        elif text[loc2] == '{':
            depth = depth + 1
        elif text[loc2] == '}':
            depth = depth - 1
        loc2 = loc2 + 1

def function_parse(text, j):
    """find the correct closing brace"""

    depth = 0
    loc2 = j
    while 1:
        if text[loc2] == '(':
            depth = depth + 1

        elif text[loc2] == ')':
            depth = depth - 1
            if depth == 0:
                return loc2
        loc2 = loc2 + 1



def get_arg_space(arg_string):
    i = arg_string.find('CellData')
    j = arg_string.find('MatData')
    k = arg_string.find('CellMatData')
    if k>=0:
        return 'CellMatData'
    if i>=0:
        return 'CellData'
    if j>=0:
        return 'MatData'

def get_arg_dim(arg_string):
    i = re.search('<([0-9])>', arg_string)
    if i is None:
        print 'error, could not determin dim: ' + arg_string
    return i.group(1)

def get_arg_neigh(arg_string):
    args = arg_parse(arg_string, 0)
    temp = {'type' : 'arg_dat',
               'type2' : 'NEIGH',
               'dim' : get_arg_dim(arg_string),
               'space' : get_arg_space(arg_string),
               'dat' : args[0],
               'acc': 'READ',
               'opt': 0,
               'arg_qtype': arg_string[0:arg_string.find('(')],
               'sten': args[1]}
    return temp

def get_arg_out(arg_string):
    args = arg_parse(arg_string, 0)
    temp = {'type' : 'arg_dat',
               'type2' : 'OUT',
               'dim' : get_arg_dim(arg_string),
               'space' : get_arg_space(arg_string),
               'dat' : args[0],
               'acc': 'WRITE',
               'opt': 0,
               'arg_qtype': arg_string[0:arg_string.find('(')],
               'sten': '0'}
    return temp

def get_arg_in(arg_string):
    args = arg_parse(arg_string, 0)
    temp = {'type' : 'arg_dat',
               'type2' : 'IN',
               'dim' : get_arg_dim(arg_string),
               'space' : get_arg_space(arg_string),
               'dat' : args[0],
               'acc': 'READ',
               'opt': 0,
               'arg_qtype': arg_string[0:arg_string.find('(')],
               'sten': '0'}
    return temp

def get_arg_reduce(arg_string):
    args = arg_parse(arg_string, 0)
    temp = {'type' : 'arg_reduce',
               'type2' : 'REDUCE',
               'dim' : get_arg_dim(arg_string),
               'space' : get_arg_space(arg_string),
               'dat' : args[1],
               'acc': args[0],
               'arg_qtype': arg_string[0:arg_string.find('(')],
               'opt': 0}
    return temp

def get_arg_gbl(arg_string):
    args = arg_parse(arg_string, 0)
    if arg_string.find('FREE_SCALAR') == 0:
        type2 = 'FREE_SCALAR'
        dim = 1
    elif arg_string.find('FREE_ARRAY') == 0:
        type2 = 'FREE_ARRAY'
        dim = args[1]

    temp = {'type' : 'arg_gbl',
               'type2' : type2,
               'dim' : dim,
               'space' : '',
               'dat' : args[0],
               'acc': 'READ',
               'arg_qtype': arg_string[0:arg_string.find('(')],
               'opt': 0}
    return temp


def get_arg_idx(arg_string):
    args = arg_parse(arg_string, 0)
    temp_idx = {'type': 'arg_idx',
                'type2': 'INDEX',
                'space' : '',
               'arg_qtype': arg_string[0:arg_string.find('(')],
                'dim' : get_arg_dim(arg_string)}
    return temp_idx

def mm_par_loop_parse(text):
  """Parsing for mm_par_loop calls"""

  loops = []

  text = comment_remover(text)
  search = ".compute"
  i = text.find(search)
  while i > -1:
      arg_list = arg_parse(text, i)
      # parse arguments in par loop
      temp_args = []
      num_args = len(arg_list)-1

      # parse each op_arg_dat
      for arg in arg_list:
        if arg.find("INDEX") == 0:
            temp_args.append(get_arg_idx(arg))
        elif arg.find("IN") == 0:
            temp_args.append(get_arg_in(arg))
        elif arg.find("OUT") == 0:
            temp_args.append(get_arg_out(arg))
        elif arg.find("NEIGH") == 0:
            temp_args.append(get_arg_neigh(arg))
        elif arg.find("FREE") == 0:
            temp_args.append(get_arg_gbl(arg))
        elif arg.find("REDUCE") == 0:
            temp_args.append(get_arg_reduce(arg))

      space = ''
      dim = 0
      for arg in temp_args:
          if arg['space'] == 'CellMatData':
              space = 'CellMatData'
          dim = max(dim, int(arg['dim']))
      if space == '':
          for arg in temp_args: #assumption: cannot have CellData and MatData and no CellMatData
              if arg['space'] == 'MatData' or arg['space'] == 'CellData':
                  space = arg['space']


      wspace1 = text[0:i].rfind(' ')
      wspace1 = max(wspace1, text[0:i].rfind('\n'))
      wspace1 = max(wspace1, text[0:i].rfind('\t'))
      computation = text[wspace1+1:i]

      if arg_list[0].find('[')>=0:
        name = 'anonymusAt'+str(i)
        lam = arg_list[0]
      else:
        name = arg_list[0]
        matches = re.findall('\\b'+arg_list[0]+'\\b\\s*=\\s*\[',text, re.MULTILINE)
        if len(matches) == 0:
            print 'Error, kernel ' + arg_list[0] + ' implementation not found'
        end = matches[-1].end(0)
        end = text.find(']',end)
        print 'TODO find LAMBDA or function'

      temp = {'loc': wspace1+1,
            'name': name,
            'computation' : computation,
            'dim' : dim,
            'lam' : lam,
            'space': space,
            'args': temp_args,
            'arg_list': arg_list[1:],
            'nargs': num_args}

      loops.append(temp)

      i = text.find(search, i + len(search))

  return (loops)

def main(source_files):

  accs_labels = ['READ', 'WRITE', 'RW', 'INC',
                    'MAX', 'MIN']
  if not source_files:
    raise ValueError("No source files specified.")

  # declare constants

  nkernels = 0
  kernels = []
  kernels_in_files = []

  #
  # loop over all input source files
  #
  
  # Find the macros defiend in the source files 
  for a in range(0, len(source_files)):
        src_file = str(source_files[a])
        f = open(src_file, 'r')
        text = f.read()
        f.close()


  kernels_in_files = [[] for _ in range(len(source_files))]
  for a in range(0, len(source_files)):
      print 'processing file ' + str(a) + ' of ' + str(len(source_files)) + \
            ' ' + str(source_files[a])

      src_file = str(source_files[a])
      f = open(src_file, 'r')
      text = f.read()

      #get rid of all comments
      text = remove_trailing_w_space(comment_remover(text))

      loops = mm_par_loop_parse(text)

      for i in range(0, len(loops)):
        name = loops[i]['name']
        nargs = loops[i]['nargs']
        dim   = loops[i]['dim']
        space = loops[i]['space']
        computation = loops[i]['computation']
        lam = loops[i]['lam']
        arg_list = loops[i]['arg_list']
        print '\nprocessing kernel ' + name + ' with ' + str(nargs) + ' arguments'
        print 'dim: '+ str(dim)
        print 'lambda: '+ lam

        #
        # process arguments
        #
        types = [''] * nargs
        types2 = [''] * nargs
        var = [''] * nargs
        stens = [''] * nargs
        accs = [''] * nargs
        dims = [''] * nargs
        spaces = [''] * nargs
        arg_qtypes = [''] * nargs
        opts = [0] * nargs

        for m in range(0, nargs):
          args = loops[i]['args'][m]
          arg_type = args['type']
          arg_qtypes[m] = args['arg_qtype']

          if arg_type == 'arg_dat':
            types[m] = args['type']
            types2[m] = args['type2']
            dims[m] = args['dim']
            spaces[m] = args['space']
            var[m] = args['dat']
            stens[m] = args['sten']
            opts[m] = args['opt']

            l = -1
            for l in range(0, len(accs_labels)):
                if args['acc'] == accs_labels[l]:
                  break

            if l == -1:
                print 'unknown access type for argument ' + str(m)
            else:
                accs[m] = l + 1

            print var[m]+' '+types2[m]+' '+spaces[m]+'<'+str(dims[m]) +'> '+str(stens[m])+' '+str(accs[m])


          if arg_type == 'arg_gbl' or arg_type == 'arg_reduce':
            types[m] = args['type']
            types2[m] = args['type2']
            dims[m] = args['dim']
            spaces[m] = args['space']
            var[m] = args['dat']
            opts[m] = args['opt']


            l = -1
            for l in range(0, len(accs_labels)):
                if args['acc'] == accs_labels[l]:
                    break
            if l == -1:
                print 'unknown access type for argument ' + str(m)
            else:
                accs[m] = l + 1

            print var[m]+' '+ str(dims[m]) +' '+str(accs[m])

          if arg_type.strip() == 'arg_idx':
            types[m] = args['type']
            types2[m] = args['type2']
            dims[m] = args['dim']
            print 'arg_idx'


        #
        # check for repeats
        #
        repeat = False
        rep1 = False
        rep2 = False
        which_file = -1
        for nk in range(0, nkernels):
          rep1 = kernels[nk]['name'] == name and \
            kernels[nk]['nargs'] == nargs and \
            kernels[nk]['dim'] == dim and \
            kernels[nk]['range'] == _range
          if rep1:
            rep2 = True
            for arg in range(0, nargs):
                rep2 = rep2 and \
                    kernels[nk]['stens'][arg] == stens[arg] and \
                    kernels[nk]['dims'][arg] == dims[arg] and \
                    kernels[nk]['types2'][arg] == typs[arg] and \
                    kernels[nk]['accs'][arg] == accs[arg]
            if rep2:
              print 'repeated kernel with compatible arguments: ' + \
                    kernels[nk]['name'],
              repeat = True
              which_file = nk
            else:
              print 'repeated kernel with incompatible arguments: ERROR'
              break


        #
        # output various diagnostics
        #

        ##
        ##todo -- not sure what will be interesting here
        ##

        #
        # store away in master list
        #

        if not repeat:
            nkernels = nkernels + 1
            temp = { 'arg_type':types,
                    'arg_type2':types2,
                    'arg_qtypes':arg_qtypes,
                    'arg_list':arg_list,
                     'name': name,
                    'nargs': nargs,
                    'dim': dim,
                    'dims': dims,
                    'stens': stens,
                    'var': var,
                    'accs': accs,
                    'spaces':spaces,
                    'opts':opts,
                    'space':space,
                    'computation':computation,
                    'lam':lam

            }
            kernels.append(temp)
            (kernels_in_files[a - 1]).append(nkernels - 1)
        else:
            append = 1
            for in_file in range(0, len(kernels_in_files[a - 1])):
                if kernels_in_files[a - 1][in_file] == which_file:
                    append = 0
            if append == 1:
                (kernels_in_files[a - 1]).append(which_file)


      #
      # output new source file
      #

      fid = open(src_file.split('.')[0] + '_mm.cpp', 'w')
      date = datetime.datetime.now()
      fid.write('//\n// auto-generated by mm.py\n//\n')

      loc_old = 0

      # read original file and locate header location
      loc_header = text.find('using namespace MM::')
      loc_header = text.find(';',loc_header)
      loc_header = text.find('\n',loc_header)+1
      loc_header = [loc_header]

      # get locations of all ops_par_loops
      n_loops = len(loops)
      loc_loops = [0] * n_loops
      for n in range(0, n_loops):
          loc_loops[n] = loops[n]['loc']

      #get locations of all kernel.h headder file declarations
      loc_kernel_headers = []

      locs = sorted(loc_header + loc_loops)

      # process header and loops
      for loc in range(0, len(locs)):
        if locs[loc] != -1:
            fid.write(text[loc_old:locs[loc] - 1])
            loc_old = locs[loc] - 1

        indent = ''
        ind = 0
        while 1:
            if text[locs[loc] - ind] == '\n':
                break
            indent = indent + ' '
            ind = ind + 1

        if (locs[loc] in loc_header) and (locs[loc] != -1):
            if len(kernels_in_files[a - 1]) > 0:
              fid.write('\n//\n// mm_par_loop declarations\n//\n')
            for k_iter in range(0, len(kernels_in_files[a - 1])):
                k = kernels_in_files[a - 1][k_iter]
                line = '\nvoid mm_par_loop_' + \
                        kernels[k]['name'] + '(std::string, Computation<'+str(kernels[k]['dim'])+'>&,\n'
                for n in range(0, kernels[k]['nargs']-1):
                    line = line + '  '+kernels[k]['arg_qtypes'][n]+',\n'
                line = line + '  '+kernels[k]['arg_qtypes'][kernels[k]['nargs']-1]+');\n'
                fid.write(line)

            fid.write('\n')
            loc_old = locs[loc] + 1
            continue


        if locs[loc] in loc_loops:
          indent = indent + ' ' * len('mm_par_loop')
          endofcall = function_parse(text, locs[loc])
          curr_loop = loc_loops.index(locs[loc])
          name = loops[curr_loop]['name']
          line = str(' mm_par_loop_' + name + '("' +
                     loops[curr_loop]['name'] + '", ' +
                     loops[curr_loop]['computation'] + ',\n' + indent)

          for arguments in range(0, loops[curr_loop]['nargs']):
              line = line + loops[curr_loop]['arg_list'][arguments]+',\n'+indent

          fid.write(line[0:-len(indent) - 2] + ');')

          loc_old = endofcall + 1
          continue

      fid.write(text[loc_old:])
      fid.close()
      f.close()



  #
  # finally, generate target-specific kernel files
  #


  mm_gen_mpi(str(source_files[0]), date, kernels)

if __name__ == '__main__':
    if len(sys.argv) > 1:
        main(source_files=sys.argv[1:]) # [1:] ignores the ops.py file itself.
    # Print usage message if no arguments given
    else:
        print __doc__
        sys.exit(1)
