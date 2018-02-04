#!/usr/bin/env python3
# -*-coding:utf-8-*-

__version__ = "v1.01.1"
__author__ = "LiuJia"

__str__ = """
********** MANUSCRIPT **********

A simple Python script used to calculate effective masses automatically.
This script is released under GPL v3.0 license.

AUTHOR: LiuJia
WEBSITE: https://github.com/liujiacode/EffectiveMass

USAGE:
# CMD or bash shell
python EffectiveMass.py -p "project_directory" [-s "xsd_file"] [-c "bs_xcd_file"] [-h] [-v]

1. This script would create a new file named "'project_name' Effective Mass.out"
   in the specified project directory, which contains the calculation results.
2. Option "-p" (or "--project") is required to specify the project directory (i.e. "project_directory")
   which contains your calculated model file ("xsd_file") and band structure file ("bs_xcd_file").
3. Options "-s" ("--xsd") and "-c" ("--xcd") are not essential
   if you don't change the project directory or file names.
   The default names of "xsd_file" and "bs_xcd_file" are
   "'project_name'.xsd" and "'project_name' Band Structure.xcd", respectively.
4. Options "-h" ("--help") and "-v" ("--version") are used to view the manuscript and
   check the version of EffectiveMass, respectively.

************* END *************
"""

import math
import os
import sys
from lxml import html
import re
import getopt


class Vector(object):
    def __init__(self, vector):
        if isinstance(vector, Vector):
            self._x, self._y, self._z = vector.values()
        elif isinstance(vector, tuple):
            if len(vector) == 3:
                self._x, self._y, self._z = (float(i) for i in vector)
            else:
                raise ValueError("Invalid input tuple: %s." % vector)
        else:
            raise ValueError("Invalid input vector: %s." % vector)

    def values(self):
        return self._x, self._y, self._z

    def __repr__(self):
        return str(self.values())

    def __getitem__(self, item):
        if isinstance(item, int) and 0 <= item <= 2:
            return self.values()[item]
        elif isinstance(item, str) and item in ('x', 'y', 'z'):
            return {'x': self._x, 'y': self._y, 'z': self._z}[item]
        else:
            raise ValueError("Invalid item: {}.".format(item))

    def __add__(self, other):
        if isinstance(other, Vector):
            return Vector(tuple(map(lambda x, y: x + y, self.values(), other.values())))
        else:
            raise ValueError("Invalid vector for add: {}.".format(other))

    def __sub__(self, other):
        if isinstance(other, Vector):
            return Vector(tuple(map(lambda x, y: x - y, self.values(), other.values())))
        else:
            raise ValueError("Invalid vector for sub: {}.".format(other))

    def __mul__(self, other):
        if isinstance(other, int) or isinstance(other, float):
            return Vector(tuple(map(lambda x: x * other, self.values())))
        else:
            raise ValueError("Invalid number for mul: {}.".format(other))

    def __truediv__(self, other):
        if isinstance(other, int) or isinstance(other, float):
            if other == 0.0:
                raise ValueError("Invalid truediv for zero.")
            else:
                return Vector(tuple(map(lambda x: x / other, self.values())))
        else:
            raise ValueError("Invalid number for truediv: {}.".format(other))

    def __floordiv__(self, other):
        if isinstance(other, int) or isinstance(other, float):
            if other == 0.0:
                raise ValueError("Invalid truediv for zero.")
            else:
                return Vector(tuple(map(lambda x: x // other, self.values())))
        else:
            raise ValueError("Invalid number for floordiv: {}.".format(other))

    def number_product(self, other):
        if isinstance(other, Vector):
            return sum(tuple(map(lambda x, y: x * y, self.values(), other.values())))
        else:
            raise ValueError("Invalid vector for number product: {}.".format(other))

    def vector_product(self, other):
        if isinstance(other, Vector):
            if self == Vector((0, 0, 0)) or other == Vector((0, 0, 0)):
                raise ValueError("Invalid vector product for zero vector.")
            else:
                return Vector((self._y * other['z'] - other['y'] * self._z,
                               self._z * other['x'] - other['z'] * self._x,
                               self._x * other['y'] - other['x'] * self._y))
        else:
            raise ValueError("Invalid vector for vector product: {}.".format(other))

    def __eq__(self, other):
        if isinstance(other, Vector):
            if self.values() == other.values():
                return True
            else:
                return False
        else:
            raise ValueError("Invalid vector for abs_equal: {}.".format(other))

    def __ne__(self, other):
        if isinstance(other, Vector):
            if self.values() == other.values():
                return False
            else:
                return True
        else:
            raise ValueError("Invalid vector for ne: {}.".format(other))

    def length(self):
        return math.sqrt(sum(map(lambda x: x ** 2, self.values())))


class EffectiveMass(object):
    def __init__(self, project):
        if os.path.isdir(project):
            self._project = project
        else:
            raise ValueError("Invalid project path, no project exist: %s." % project)

        self._project_name = os.path.basename(self._project)
        self._name = ' '.join(self._project_name.split(' ')[:-2])
        self._code = self._project_name.split(' ')[-2]
        self._method = self._project_name.split(' ')[-1]
        if self._code not in ["DMol3", "CASTEP"] or self._method not in ["GeomOpt", "Energy"]:
            raise ValueError("Invalid code or method of the project: %s." % self._project_name)

        output_file = self._project + os.sep + self._name + " Effective Mass"
        file_index = 1
        if os.path.isfile(output_file + ".out"):
            while os.path.isfile(output_file + '(%s).out' % str(file_index)):
                file_index += 1
            output_file = output_file + '(%s)' % str(file_index)
        self._output_file = output_file + ".out"
        with open(self._output_file, 'w') as o_file:
            o_file.write("Electron and hole effective mass calculated by %s.\n" % sys.path[0])
            o_file.write("EffectiveMass.py Version: %s.\n" % __version__)
            o_file.write("Project: %s.\n" % self._project)
            o_file.write("Code: %s\n" % self._code)
            o_file.write("Method: %s\n" % self._method)

    def output(self, string):
        with open(self._output_file, 'a') as o_file:
            o_file.write(str(string) + '\n')

    def load(self, xsd_file="default", bs_xcd_file="default"):
        if xsd_file == "default":
            xsd_file = self._project + os.sep + self._name + ".xsd"
        else:
            xsd_file = self._project + os.sep + str(xsd_file)
        if os.path.isfile(xsd_file):
            xsd_data = html.fromstring(open(xsd_file, 'r').read().encode("utf-8"))
        else:
            raise ValueError("xsd file file is not existed: %s." % xsd_file)

        if bs_xcd_file == "default":
            bs_xcd_file = self._project + os.sep + self._name + " Band Structure.xcd"
        else:
            bs_xcd_file = self._project + os.sep + str(bs_xcd_file)
        if os.path.isfile(bs_xcd_file):
            bs_xcd_data = html.fromstring(open(bs_xcd_file, 'r').read().encode("utf-8"))
        else:
            raise ValueError("band structure xsd file is not existed: %s." % bs_xcd_file)

        # Loading real lattice.
        space_group = xsd_data.xpath("//spacegroup")
        if len(space_group) != 1:
            raise ValueError("Too many(little) 'SpaceGroup' patterns: %s" % xsd_file)
        a_vector = Vector(tuple(float(i) for i in space_group[0].xpath("./@avector")[0].split(',')))
        b_vector = Vector(tuple(float(i) for i in space_group[0].xpath("./@bvector")[0].split(',')))
        c_vector = Vector(tuple(float(i) for i in space_group[0].xpath("./@cvector")[0].split(',')))
        cell_vectors = (a_vector, b_vector, c_vector)

        # Loading k points and k lines from xsd_file.
        k_path = xsd_data.xpath("//kpath")
        if len(k_path) != 1:
            raise ValueError("Too many(little) 'KPath' patterns: %s." % xsd_file)
        # k_path_children = list(int(i) for i in k_path[0].xpath("./@children")[0].split(','))
        k_path_objects = list(int(i) for i in k_path[0].xpath("./@objects")[0].split(','))

        k_points = {}
        k_names = {}
        for k_point in xsd_data.xpath("//kpoint"):
            k_point_id = int(k_point.xpath("./@id")[0])
            k_point_name = k_point.xpath("./@name")[0]
            k_point_HKL = Vector(tuple(float(i) for i in k_point.xpath("./@hkl")[0].split(',')))
            k_points[k_point_name] = k_point_HKL
            k_names[k_point_id] = k_point_name
        if not k_points:
            raise ValueError("Loading k points data failed.")

        k_lines = []
        for object_id in k_path_objects:
            k_line = xsd_data.xpath("//kline[@id=%s]" % object_id)
            if len(k_line) == 1:
                k_line_connects = tuple(int(i) for i in k_line[0].xpath("./@connects")[0].split(','))
                k_lines.append(tuple(k_names[i] for i in k_line_connects))
            else:
                raise ValueError("K line is not found. %s" % object_id)

        k_lines = tuple(k_lines)
        if not k_lines:
            raise ValueError("Loading k lines data failed.")

        # Loading band structures.
        series_2d = bs_xcd_data.xpath("//series_2d")
        if len(series_2d) != 1:
            raise ValueError("Too many(little) 'SERIES_2D' patterns: %s." % bs_xcd_file)

        numpoints = int(series_2d[0].xpath("@numpoints")[0])
        point_2d = series_2d[0].xpath("./point_2d")
        if len(point_2d) != numpoints:
            raise ValueError("Total number of 'POINT_2D' is not equal to 'NumPoints': %s." % str(numpoints))

        # valence band and conduction band data share one x axis.
        valence_band = tuple()
        conduction_band = tuple()
        x_axis = []
        next_conduction_band = False
        for each_point in point_2d:
            x, y = [float(i) for i in each_point.xpath("@xy")[0].split(',')]
            x_axis.append(x)
            if x == 1:
                break
        x_axis_len = len(x_axis)
        # Each group of band data has one waste record like "<POINT_2D XY="0,1e+308"/>".
        # So it has to add 1 when splitting band data into groups.
        rx_axis_len = x_axis_len + 1
        if numpoints % rx_axis_len != 0:
            raise ValueError("Invalid band point groups.")

        for band_index in range(int(numpoints / rx_axis_len)):
            y_axis = []
            for point_index in range(x_axis_len):
                x, y = [float(i) for i in
                        point_2d[rx_axis_len * band_index + point_index].xpath("@xy")[0].split(',')]
                if x == x_axis[point_index]:
                    y_axis.append(y)
                else:
                    raise ValueError("Invalid x of point %s." % (str(rx_axis_len * band_index + point_index)))

            dis_x, dis_y = [float(i) for i in
                            point_2d[rx_axis_len * band_index + x_axis_len].xpath("@xy")[0].split(',')]
            if dis_x != 0:
                raise ValueError(
                    "Invalid dis_x of last discarded point %s." % str(rx_axis_len * band_index + x_axis_len))

            if next_conduction_band:
                conduction_band = tuple(y_axis)
                break
            if max(y_axis) == 0:
                valence_band = tuple(y_axis)
                next_conduction_band = True

        x_axis = tuple(x_axis)
        if not (valence_band and conduction_band):
            raise ValueError("Loading valence band and conduction band data failed.")

        # Loading other paras.
        cart_2d_plot = bs_xcd_data.xpath("//cartesian_2d_plot")
        if len(cart_2d_plot) != 1:
            raise ValueError("Too many(little) 'CARTESIAN_2D_PLOT' patterns: %s." % bs_xcd_file)
        plot_title = cart_2d_plot[0].xpath("./title/@text")[0]
        title_re = re.compile(r"^%s Band Structure\\n Band gap is (\d+\.\d+) (Ha|eV)$" % self._code)
        title_match = title_re.match(plot_title)
        if title_match:
            band_gap_num, band_gap_unit = title_match.groups()
            band_gap_num = float(band_gap_num)
        else:
            raise ValueError("Unmatched title: %s." % plot_title)
        calculated_band_gap_num = min(conduction_band) - max(valence_band)
        if abs(calculated_band_gap_num - band_gap_num) > 0.001:
            raise ValueError("Invalid calculated band gap: %s." % calculated_band_gap_num)

        # Loading k path labels (high symmetry points in band structure files).
        num_labels = int(cart_2d_plot[0].xpath("./cartesian_axes_2d/axis_tick_labels/@numlabels")[0])
        axis_tick_labels = cart_2d_plot[0].xpath("./cartesian_axes_2d/axis_tick_labels/axis_tick_label")
        if len(axis_tick_labels) != num_labels:
            raise ValueError("Too many(little) axis tick labels: %s." % num_labels)
        k_path_labels = []
        for axis_tick_label in axis_tick_labels:
            position = float(axis_tick_label.xpath("./@position")[0])
            text = axis_tick_label.xpath("./@text")[0]
            if text in k_points:
                k_path_labels.append((text, position))
            else:
                raise ValueError("Point %s: %s is not found in k points list." % (text, position))

        k_path_labels = tuple(k_path_labels)
        if not k_path_labels:
            raise ValueError("Loading k path labels data failed.")

        # The sequence of k_lines should match with k_path_labels.
        for index in range(len(k_lines)):
            start = k_lines[index][0]
            end = k_lines[index][1]
            if start != k_path_labels[index][0] or end != k_path_labels[index + 1][0]:
                raise ValueError("Unmatched k lines and k path labels.")

        self.output('')
        self.output("*" * 20 + "Loading Start" + "*" * 20)
        self.output("Cell Vectors: %s." % str(cell_vectors))
        self.output("K Points: %s." % str(k_points))
        self.output("K Lines: %s." % str(k_lines))
        self.output("K Path Labels: %s." % str(k_path_labels))
        self.output("Band gap: %s (%s)." % (band_gap_num, band_gap_unit))
        self.output('')
        self.output("X Axis:\n    %s" % '\n    '.join((str(i) for i in x_axis)))
        self.output('')
        self.output("Valence Band (%s):\n    %s" % (band_gap_unit, '\n    '.join((str(i) for i in valence_band))))
        self.output('')
        self.output("Conduction Band (%s):\n    %s" % (band_gap_unit, '\n    '.join((str(i) for i in conduction_band))))
        self.output("*" * 20 + "Loading End" + "*" * 20 + '\n')

        k_path_paras = (k_points, k_lines, k_path_labels)
        band_paras = (band_gap_num, band_gap_unit, x_axis, valence_band, conduction_band)
        return cell_vectors, k_path_paras, band_paras

    def calculate(self, cell_vectors, k_path_paras, band_paras):
        k_points, k_lines, k_path_labels = k_path_paras
        band_gap_num, band_gap_unit, x_axis, valence_band, conduction_band = band_paras
        bohr_constant = 0.52917721092
        ha_constant = 27.2114
        # Calculate reciprocal vectors with bohr unit.
        bohr_cell_vectors = tuple(x / bohr_constant for x in cell_vectors)
        volume = (bohr_cell_vectors[0].vector_product(bohr_cell_vectors[1]).number_product(bohr_cell_vectors[2]))
        bohr_reciprocal_vectors = (
            bohr_cell_vectors[1].vector_product(bohr_cell_vectors[2]) * (2 * math.pi / volume),
            bohr_cell_vectors[2].vector_product(bohr_cell_vectors[0]) * (2 * math.pi / volume),
            bohr_cell_vectors[0].vector_product(bohr_cell_vectors[1]) * (2 * math.pi / volume),
        )

        # Calculate abs k points under reciprocal lattice (bohr_reciprocal_vectors).
        cart_k_points = {}
        for index in k_points.keys():
            cart_vector = Vector((0, 0, 0))
            for i in range(3):
                cart_vector = cart_vector + bohr_reciprocal_vectors[i] * k_points[index][i]
            cart_k_points[index] = cart_vector

        # Calculate the abs length of each k line.
        k_line_lengths = []
        for k_line in k_lines:
            start = cart_k_points[k_line[0]]
            end = cart_k_points[k_line[1]]
            k_line_lengths.append(Vector(end - start).length())

        # Calculate the abs labels of each k point along k lines.
        cart_k_path_labels = [k_path_labels[0]]
        for index in range(len(k_line_lengths)):
            abs_length = sum(k_line_lengths[:index + 1])
            cart_k_path_labels.append((k_path_labels[index + 1][0], abs_length))

        # Get the relative position of each k point in axis.
        kp_in_x_axis = []
        for k_path_label in k_path_labels:
            label = k_path_label[1]
            position = -1
            previous = position
            for x_index in range(len(x_axis)):
                if abs(label - x_axis[x_index]) < 0.0001:
                    position = x_index
                    break
            if position >= 0 and position > previous:
                kp_in_x_axis.append((k_path_label[0], position))
            else:
                raise ValueError("Cannot find the position of k path label: %s." % str(k_path_label))
        if len(kp_in_x_axis) != len(k_path_labels):
            raise ValueError("Too many(little) k points in x axis: %s." % str(kp_in_x_axis))

        # Calculate the absolute x axis values.
        cart_x_axis = [-1] * len(x_axis)
        for index in range(len(kp_in_x_axis) - 1):
            start_label = kp_in_x_axis[index][0]
            end_label = kp_in_x_axis[index + 1][0]
            start_index = kp_in_x_axis[index][1]
            end_index = kp_in_x_axis[index + 1][1]
            start_value = x_axis[start_index]
            end_value = x_axis[end_index]
            if start_label == cart_k_path_labels[index][0]:
                cart_start_value = cart_k_path_labels[index][1]
            else:
                raise ValueError("Unmatched start label: %s." % start_label)
            if end_label == cart_k_path_labels[index + 1][0]:
                cart_end_value = cart_k_path_labels[index + 1][1]
            else:
                raise ValueError("Unmatched end label: %s." % end_label)
            rate = (cart_end_value - cart_start_value) / (end_value - start_value)
            for x_index in range(end_index - start_index + 1):
                abs_x_index = start_index + x_index
                abs_x_value = cart_start_value + (x_axis[abs_x_index] - start_value) * rate
                if cart_x_axis[abs_x_index] == -1:
                    cart_x_axis[abs_x_index] = abs_x_value
                else:
                    if abs(cart_x_axis[abs_x_index] - abs_x_value) > 0.0001:
                        raise ValueError("Conversion failed (1).")

        for value in cart_x_axis:
            if value < 0:
                raise ValueError("Conversion failed (2).")

        for index in range(len(cart_x_axis) - 1):
            if cart_x_axis[index] >= cart_x_axis[index + 1]:
                raise ValueError("Conversion failed (3).")

        self.output('')
        self.output("*" * 20 + "Calculating Start" + "*" * 20)
        # self.output("Bohr Constant: %s (atom/bohr)." % str(bohr_constant))
        self.output("Volume: %s (bohr^3)." % volume)
        self.output("Bohr Reciprocal Vectors %s (bohr^-1)." % str(bohr_reciprocal_vectors))
        # self.output("Cartesian K Points: %s." % str(cart_k_points))
        # self.output("Cartesian K Line Length: %s." % str(k_line_lengths))
        # self.output("Cartesian K Path Labels: %s." % str(cart_k_path_labels))
        # self.output("Position of Each K Point in X Axis: %s." % str(kp_in_x_axis))
        self.output("Cartesian X Axis (bohr^-1):\n    %s" % '\n    '.join(str(i) for i in cart_x_axis))
        self.output("*" * 20 + "Calculating End" + "*" * 20)

        # This script can only calculate simple band structures.
        if kp_in_x_axis[0][0] != kp_in_x_axis[-1][0]:
            self.output('')
            self.output("Error: The band structure is too complex for this calculation.")
            self.output("    The start point is not coincident with the end point.")
            raise ValueError("Too complex band structure.")
        pos_list = []
        for pos in kp_in_x_axis[:-1]:
            if pos[0] not in pos_list:
                pos_list.append(pos[0])
            else:
                self.output('')
                self.output("Error: The band structure is too complex for this calculation.")
                self.output(
                    "    Point %s (not the start or end point) is shown up twice or more in the band structure." % pos[
                        0])
                raise ValueError("Too complex band structure.")

        # Effective mass of electron on valence band.
        e_value = max(valence_band)
        e_position = valence_band.index(e_value)
        if e_value != 0.0:
            raise ValueError("The band structure is broken. Please re-calculate the band structure.")

        # Two directions if e_position is a k point.
        if e_position in tuple(i[1] for i in kp_in_x_axis):
            e_index = tuple(i[1] for i in kp_in_x_axis).index(e_position)
            # First direction.
            if e_index == 0:
                direction_1 = (kp_in_x_axis[-2], kp_in_x_axis[-1], len(cart_x_axis) - 1)
            else:
                direction_1 = (kp_in_x_axis[e_index - 1], kp_in_x_axis[e_index], e_position)
            # Second direction.
            if e_index == len(kp_in_x_axis) - 1:
                direction_2 = (kp_in_x_axis[0], kp_in_x_axis[1], 0)
            else:
                direction_2 = (kp_in_x_axis[e_index], kp_in_x_axis[e_index + 1], e_position)
            e_directions = (direction_1, direction_2)
        else:
            e_directions = tuple()
            for e_index in range(len(kp_in_x_axis) - 1):
                if kp_in_x_axis[e_index][1] < e_position < kp_in_x_axis[e_index + 1][1]:
                    e_directions = ((kp_in_x_axis[e_index], kp_in_x_axis[e_index + 1], e_position),)
                    break
            if not e_directions:
                raise ValueError("Cannot find the electron directions.")

        e_effective_masses = []
        for e_direction in e_directions:
            if e_direction[1][1] - e_direction[0][1] < 4:
                raise ValueError("Too few points between direction: %s." % str(e_direction))
            dir_name = "%s%s" % (e_direction[0][0], e_direction[1][0])
            start_pos = e_direction[0][1]
            end_pos = e_direction[1][1]
            point_pos = e_direction[2]
            if point_pos == start_pos:
                p_1 = start_pos
                p_2 = start_pos + 1
                p_3 = start_pos + 2
            elif point_pos == end_pos:
                p_1 = end_pos - 2
                p_2 = end_pos - 1
                p_3 = end_pos
            else:
                p_1 = point_pos - 1
                p_2 = point_pos
                p_3 = point_pos + 1
            x_1 = cart_x_axis[p_1]
            x_2 = cart_x_axis[p_2]
            x_3 = cart_x_axis[p_3]
            y_1 = valence_band[p_1]
            y_2 = valence_band[p_2]
            y_3 = valence_band[p_3]
            valence_band_unit = band_gap_unit
            if valence_band_unit == 'eV':
                y_1 = y_1 / ha_constant
                y_2 = y_2 / ha_constant
                y_3 = y_3 / ha_constant
                valence_band_unit = 'Ha'
            if valence_band_unit != 'Ha':
                raise ValueError("Invalid valence band unit: %s." % valence_band_unit)

            d_2 = (y_3 - y_2) / (x_3 - x_2)
            d_1 = (y_2 - y_1) / (x_2 - x_1)
            dd = (d_2 - d_1) / ((x_3 - x_1) / 2)
            mass = 1 / dd
            e_effective_masses.append((dir_name, mass))

        if not e_effective_masses:
            raise ValueError("Calculating electron effective mass failed.")

        self.output('')
        self.output("*" * 20 + "Results Start" + "*" * 20)
        self.output("+++ me is the effective mass of free electrons. +++")
        self.output('')
        self.output('Valence Max value: %s %s.' % (e_value, band_gap_unit))
        self.output("Mass Position: %s." % e_position)
        self.output('')
        dir_content = ""
        for e_dir in e_directions:
            dir_content += "    %s(%s) -> %s(%s)\n" % (e_dir[0][0], e_dir[0][1], e_dir[1][0], e_dir[1][1])
        self.output("Electron Directions:\n%s" % dir_content)
        mass_content = ""
        for e_mass in e_effective_masses:
            mass_content += "    Direction: %s, Mass: %s me\n" % (e_mass[0], e_mass[1])
        self.output("Electron Effective Mass:\n%s" % mass_content)

        # Effective mass of hole on conduction band.
        h_value = min(conduction_band)
        h_position = conduction_band.index(h_value)

        # Two directions if h_position is a k point.
        if h_position in tuple(i[1] for i in kp_in_x_axis):
            h_index = tuple(i[1] for i in kp_in_x_axis).index(h_position)
            # First direction.
            if h_index == 0:
                direction_1 = (kp_in_x_axis[-2], kp_in_x_axis[-1], len(cart_x_axis) - 1)
            else:
                direction_1 = (kp_in_x_axis[h_index - 1], kp_in_x_axis[h_index], h_position)
            # Second direction.
            if h_index == len(kp_in_x_axis) - 1:
                direction_2 = (kp_in_x_axis[0], kp_in_x_axis[1], 0)
            else:
                direction_2 = (kp_in_x_axis[h_index], kp_in_x_axis[h_index + 1], h_position)
            h_directions = (direction_1, direction_2)
        else:
            h_directions = tuple()
            for h_index in range(len(kp_in_x_axis) - 1):
                if kp_in_x_axis[h_index][1] < h_position < kp_in_x_axis[h_index + 1][1]:
                    h_directions = ((kp_in_x_axis[h_index], kp_in_x_axis[h_index + 1], h_position),)
                    break
            if not h_directions:
                raise ValueError("Cannot find the electron directions.")

        h_effective_masses = []
        for h_direction in h_directions:
            if h_direction[1][1] - h_direction[0][1] < 4:
                raise ValueError("Too few points between direction: %s." % str(h_direction))
            dir_name = "%s%s" % (h_direction[0][0], h_direction[1][0])
            start_pos = h_direction[0][1]
            end_pos = h_direction[1][1]
            point_pos = h_direction[2]
            if point_pos == start_pos:
                p_1 = start_pos
                p_2 = start_pos + 1
                p_3 = start_pos + 2
            elif point_pos == end_pos:
                p_1 = end_pos - 2
                p_2 = end_pos - 1
                p_3 = end_pos
            else:
                p_1 = point_pos - 1
                p_2 = point_pos
                p_3 = point_pos + 1
            x_1 = cart_x_axis[p_1]
            x_2 = cart_x_axis[p_2]
            x_3 = cart_x_axis[p_3]
            y_1 = conduction_band[p_1]
            y_2 = conduction_band[p_2]
            y_3 = conduction_band[p_3]
            conduction_band_unit = band_gap_unit
            if conduction_band_unit == 'eV':
                y_1 = y_1 / ha_constant
                y_2 = y_2 / ha_constant
                y_3 = y_3 / ha_constant
                conduction_band_unit = 'Ha'
            if conduction_band_unit != 'Ha':
                raise ValueError("Invalid conduction band unit: %s." % conduction_band_unit)

            d_2 = (y_3 - y_2) / (x_3 - x_2)
            d_1 = (y_2 - y_1) / (x_2 - x_1)
            dd = (d_2 - d_1) / ((x_3 - x_1) / 2)
            mass = 1 / dd
            h_effective_masses.append((dir_name, mass))

        if not h_effective_masses:
            raise ValueError("Calculating hole effective mass failed.")

        self.output('')
        self.output('Conduction Min value: %s %s.' % (h_value, band_gap_unit))
        self.output("Mass Position: %s." % h_position)
        self.output('')
        dir_content = ""
        for h_dir in h_directions:
            dir_content += "    %s(%s) -> %s(%s)\n" % (h_dir[0][0], h_dir[0][1], h_dir[1][0], h_dir[1][1])
        self.output("Hole Directions:\n%s" % dir_content)
        mass_content = ""
        for h_mass in h_effective_masses:
            mass_content += "    Direction: %s, Mass: %s me\n" % (h_mass[0], h_mass[1])
        self.output("Hole Effective Mass:\n%s" % mass_content)
        self.output("*" * 20 + "Results End" + "*" * 20)


if __name__ == '__main__':
    opts, args = getopt.getopt(sys.argv[1:], "p:s:c:hv", ["project=", "xsd=", "bs_xcd=", "help", "version"])
    project_directory = ""
    model_xsd_file = "default"
    bandstr_xcd_file = "default"

    for op, value in opts:
        if op in ("-h", "--help"):
            print(__str__)
            exit()
        elif op in ("-v", "--version"):
            print("EffectiveMass %s" % __version__)
            exit()
        elif op in ("-p", "--project"):
            project_directory = value
        elif op in ("-s", "--xsd"):
            model_xsd_file = value
        elif op in ("-c", "--bs_xcd"):
            bandstr_xcd_file = value
        else:
            raise ValueError("Invalid option: %s." % op)

    if not project_directory:
        raise ValueError("Project directory is empty.")

    em = EffectiveMass(project_directory)
    c_paras, k_paras, b_paras = em.load(xsd_file=model_xsd_file, bs_xcd_file=bandstr_xcd_file)
    em.calculate(cell_vectors=c_paras, k_path_paras=k_paras, band_paras=b_paras)
