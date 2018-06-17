#!/usr/bin/env python3
# -*-coding:utf-8-*-

__version__ = "v1.0.0"
__author__ = "Jia Liu"

__str__ = """
********** HELP MANUSCRIPT **********

[INTRODUCTION]
  Calculate electron density difference manually. (for DMol3)
  For CASTEP please use calculation property "Electron density difference" with "Sets of atoms".
  Note that the total density and densities of different sets should be calculated before running this script.
  This script is released under GPL v3.0 license.
  Current version: %s
[AUTHOR] %s
[WEBSITE] https://github.com/liujiacode/EffectiveMass

[USAGE]
  >> python DensitySub.py -t "total_file" -s "sub_files" -a "add_files" -o "output_file" [-h] [-v]

[PARAMETERS]
total electron density file ("sub_files" and "add_files"):
  Total electron density file.
  Note that the density files are usually hidden files with ".grd" extensions.

electron density files ("sub_files" and "add_files"):
  Electron density files for subtracting and adding from/to total density.
  The file paths should be separated by ';' without spaces.

output density file ("output_file")
  For output.

**************** END ****************
""" % (__version__, __author__)

import os
import sys
import getopt


class DensitySub(object):
    def __init__(self, total, subs, adds):
        if os.path.isfile(total):
            self._total = open(total,'r')
        else:
            raise ValueError("Invalid total density file: %s." % total)
        
        if isinstance(subs, list) or isinstance(subs, tuple):
            self._subs = []
            for sub in subs:
                if os.path.isfile(sub):
                    self._subs.append(open(sub, 'r'))
                else:
                    raise ValueError("Invalid sub density file: %s." % sub)
        
        if isinstance(adds, list) or isinstance(adds, tuple):
            self._adds = []
            for add in adds:
                if os.path.isfile(add):
                    self._adds.append(open(add, 'r'))
                else:
                    raise ValueError("Invalid add density file: %s." % add)

    def calculate(self):
        total_data = self._total.read().split('\n')
        while '' in total_data:
            total_data.remove('')

        for sub in self._subs:
            sub_data = sub.read().split('\n')
            while '' in sub_data:
                sub_data.remove('')
            
            if sub_data[:5] != total_data[:5]:
                raise ValueError("The crystal infos for sub are mismatching.")
            if len(sub_data) != len(total_data):
                raise ValueError("The points length for sub is mismatching.")
            
            for poi in range(len(total_data) - 5):
                total_data[5 + poi] = float(total_data[5 + poi]) - float(sub_data[5 + poi])
        
        for add in self._adds:
            add_data = add.read().split('\n')
            while '' in add_data:
                add_data.remove('')
            
            if add_data[:5] != total_data[:5]:
                raise ValueError("The crystal infos for add are mismatching.")
            if len(add_data) != len(total_data):
                raise ValueError("The points length for add is mismatching.")
            
            for poi in range(len(total_data) - 5):
                total_data[5 + poi] = float(total_data[5 + poi]) + float(add_data[5 + poi])
        
        return total_data
        
    def write(self, output_path):
        total_data = self.calculate()
        with open(output_path, 'w') as output:
            for data in total_data:
                output.write(str(data) + '\n')

                
if __name__ == '__main__':
    opts, args = getopt.getopt(sys.argv[1:], "t:s:a:o:hv", ["total=", "subs=", "adds=", "output", "help", "version"])
    total_file = ""
    sub_files = []
    add_files = []
    output_file = ""

    for op, value in opts:
        if op in ("-h", "--help"):
            print(__str__)
            exit()
        elif op in ("-v", "--version"):
            print("DensitySub %s" % __version__)
            exit()
        elif op in ("-t", "--total"):
            total_file = value
        elif op in ("-s", "--subs"):
            sub_files = value.split(';')
        elif op in ("-a", "--adds"):
            add_files = value.split(';')
        elif op in ("-o", "--output"):
            output_file = value
        else:
            raise ValueError("Invalid option: %s." % op)

    ds = DensitySub(total=total_file, subs=sub_files, adds=add_files)
    ds.write(output_file)
