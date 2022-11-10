#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import re
from sys import argv, exit, stdout, stderr
from xml.etree.ElementTree import Element, SubElement, Comment, tostring
from xml.dom import minidom
import numpy as np

def prettify(elem):
    """Return a pretty-printed XML string for the Element.
    """
    rough_string = tostring(elem, 'utf-8')
    reparsed = minidom.parseString(rough_string)
    return reparsed.toprettyxml(indent="  ")

def five_columns(data, format="%0.4e "):
    s = ""
    for i in range(len(data)//5):
      s += format % data[5*i]
      s += format % data[5*i+1]
      s += format % data[5*i+2]
      s += format % data[5*i+3]
      s += (format % data[5*i+4])
      if (5*i+4 != len(data)-1):
        s += "\n"
#
    if (len(data)%5 != 0):
      for i in range(5*(len(data)//5),len(data),1):
        s += format % data[i]
    return s

def read_card(nsysi, names, values, types):
    """
    Reads cards from a file.
    """
    if (types[0] != 't'):
        l = nsysi.readline()
        l = l.split('\n')[0]
        line = l.split('/')[0]
        pattern = re.compile('\s+')
        line = re.sub(pattern, ' ', line)
        pattern = re.compile('^\s')
        line = re.sub(pattern, '', line)
        pattern = re.compile('\s+$')
        line = re.sub(pattern, '', line)
        line = line.split(' ')
        d = dict(zip(names, values))
        for i in range(min(len(names),len(line))):
            if types[i] == 'i':
                try:
                    d[names[i]] = int(line[i])
                except (TypeError, ValueError):
                    stderr.write("".join(['Error while loading "', names[i], '".\n']))
                    stderr.write("".join(['Integer expected and [',line[i],'] cannot be converted.\n']))
                    stderr.write("".join(['Line: "', l,'".\n']))
                    exit(0)
            elif types[i] == 'r':
                try:
                    d[names[i]] = float(line[i])
                except (TypeError, ValueError):
                    stderr.write("".join(['Error while loading "', names[i], '".\n']))
                    stderr.write("".join(['Real expected and [',line[i],'] cannot be converted.\n']))
                    stderr.write("".join(['Line: "', l,'".\n']))
                    exit(0)
    else:
        l = nsysi.readline()
        l = l.split('\n')[0]
        line = l.split('/')[0]
        d = dict(zip(names, values))
        d[names[0]] = line
    return d
    
def read_floatarray(nsysi, nvalues, name):
    """
    Reads a float array from a file.
    """
    list = []
    i = 0
    while i < nvalues:
        l = nsysi.readline()
        l = l.split('\n')[0]
        line = l.split('/')[0]
        pattern = re.compile('\s+')
        line = re.sub(pattern, ' ', line)
        pattern = re.compile('^\s')
        line = re.sub(pattern, '', line)
        pattern = re.compile('\s+$')
        line = re.sub(pattern, '', line)
        line = line.split(' ')
        for value in line:
            if value != '':
                try:
                    list.append(float(value))
                    i = i + 1
                except (TypeError, ValueError):
                    stderr.write("".join(['Error while loading element ', str(i), ' of array "',name,'".\n']))
                    stderr.write("".join(['Real expected and [',value,'] cannot be converted.\n']))
                    stderr.write("".join(['Line: "', l,'".\n']))
                    exit(0)
    return list

def leapr2xml(leapr_input, xml_output=""):
    """
    Reads a LEAPR input and creates an XML file. If xml_output is given, the file is saved to disk.
    If not, the input is returned as a string.
    """

#
# Load default values
#
    card_names =  [['nout'],
                   ['title'],
                   ['ntempr','iprint','nphon'],
                   ['mat','za','isabt','ilog'],
                   ['awr','spr','npr','iel','ncold','nsk'],
                   ['nss', 'b7','aws', 'sps','mss'],
                   ['nalpha', 'nbeta', 'lat']]
    card_values = [[0],
                   [""],
                   [0,1,100],
                   [0,0,0,0],
                   [0,0,0,0,0,0],
                   [0,0,0,0,0],
                   [0,0,0]]
    card_type = [['i'],
                   ['t'],
                   ['i','i','i'],
                   ['i','r','i','i'],
                   ['r','r','i','i','i','i'],
                   ['i','i','r','r','i'],
                   ['i','i','i']]
    tcard_names =  [['temp'],
                    ['delta', 'ni'],
                    ['twt','c','tbeta'],
                    ['nd'],
                    ['nka', 'dka'],
                    ['cfrac']]
    tcard_values = [[0],
                    [0,0],
                    [0,0,0],
                    [0],
                    [0, 0],
                    [0]]
    tcard_type =   [['r'],
                    ['r','i'],
                    ['r','r','r'],
                    ['i'],
                    ['i', 'r'],
                    ['r']]
        
    try:
        nsysi = open(leapr_input)
    except IOError:
        print("Error: ", leapr_input, " does not exist")
        raise SystemExit
    line = nsysi.readline()
    line = line.split('\n')[0]
    line2 = line.upper()
    if (line2.find("LEAPR") == -1):
        print("Error: ", leapr_input, " is not a LEAPR input")
        print("LEAPR expected in the first line but the file starts with")
        print(line)
        raise SystemExit

    params = {}
    for i in range(len(card_names)):
        d = read_card(nsysi, card_names[i], card_values[i], card_type[i])
        params.update(d)
    params['alphas'] = read_floatarray(nsysi, params['nalpha'], 'alphas')
    params['betas']  = read_floatarray(nsysi, params['nbeta'], 'betas')
    tparams = []
    for t in range(params['ntempr']):
        tparams.append({})
        d = read_card(nsysi, tcard_names[0], tcard_values[0], tcard_type[0])
        tparams[t].update(d)
        if tparams[t]['temp'] > 0:
            d = read_card(nsysi, tcard_names[1], tcard_values[1], tcard_type[1])
            tparams[t].update(d)
            tparams[t]['rho'] = read_floatarray(nsysi, tparams[t]['ni'], 'rho')
            d = read_card(nsysi, tcard_names[2], tcard_values[2], tcard_type[2])
            tparams[t].update(d)
            d = read_card(nsysi, tcard_names[3], tcard_values[3], tcard_type[3])
            tparams[t].update(d)
            if tparams[t]['nd'] > 0:
                tparams[t]['bdel'] = read_floatarray(nsysi, tparams[t]['nd'], 'bdel')
                tparams[t]['adel'] = read_floatarray(nsysi, tparams[t]['nd'], 'adel')
            if params['ncold'] > 0 or params['nsk'] > 0:
                d = read_card(nsysi, tcard_names[4], tcard_values[4], tcard_type[4])
                tparams[t].update(d)
                tparams[t]['skappa'] = read_floatarray(nsysi, tparams[t]['nka'], 'skappa')
            if ((params['ncold'] == 0) and (params['nsk'] > 0)):
            #
            # cfrac is hard coded for ncold
            #
                d = read_card(nsysi, tcard_names[5], tcard_values[5], tcard_type[5])
                tparams[t].update(d)
        else:
            tparams[t]['temp'] = abs(tparams[t]['temp'])
            tparams[t]['delta'] = tparams[t-1]['delta'] 
            tparams[t]['ni'] = tparams[t-1]['ni'] 
            tparams[t]['rho'] = tparams[t-1]['rho'] 
            tparams[t]['twt'] = tparams[t-1]['twt'] 
            tparams[t]['c'] = tparams[t-1]['c'] 
            tparams[t]['tbeta'] = tparams[t-1]['tbeta'] 
            tparams[t]['nd'] = tparams[t-1]['nd'] 
            if tparams[t]['nd'] > 0:
                tparams[t]['bdel'] = tparams[t-1]['bdel']
                tparams[t]['adel'] = tparams[t-1]['adel']
            if params['ncold'] > 0 or params['nsk'] > 0:
                tparams[t]['nsk'] = tparams[t-1]['nsk']
                tparams[t]['dsk'] = tparams[t-1]['dsk']
                tparams[t]['skappa'] = tparams[t-1]['skappa']
            if ((params['ncold'] == 0) and (params['nsk'] > 0)):
                tparams[t]['cfrac'] = tparams[t-1]['cfrac']
    comments = ""
    while True:
        l = nsysi.readline()
        if (not l):
            break
        l2 = l.split('\n')[0]
        l2 = l2.split('/')[0]
        if ((l2=="") or (l2=="stop")):
            break
        comments = comments + l
    
    model = Element('model')
        
    leapr = SubElement(model, 'leapr')
    
    for key in params:
        elem = SubElement(leapr, key)
        if type(params[key]) == int:
            elem.text = "%i" % (params[key])
        if type(params[key]) == float:
            elem.text = "%0.4e" % (params[key])
        if type(params[key]) == str:
            elem.text = "%s" % (params[key])
        if key in ['alphas', 'betas']:
#            s = np.array2string(np.array(params[key]),  max_line_width=55, formatter={'float_kind':lambda x: "%.4e" % x})
#            s = re.sub('[\[\]]', '', s)
#            s = re.sub('\n ', ' \n', s)
            s = five_columns(np.array(params[key]))
            elem.text = s
    
    for t in range(params['ntempr']):
        temperature = SubElement(leapr, 'temperature')
        for key in tparams[t]:
            elem = SubElement(temperature, key)
            if key in ['rho', 'adel', 'bdel', 'skappa']:
#                    s = np.array2string(np.array(tparams[t][key]),  max_line_width=55, formatter={'float_kind':lambda x: "%.4e" % x})
#                    s = re.sub('[\[\]]', '', s)
#                    s = re.sub('\n ', ' \n', s)
                    s = five_columns(np.array(tparams[t][key]))
                    elem.text = s
            else:
                if type(tparams[t][key]) == int:
                    elem.text = "%i" % (tparams[t][key])
                if type(tparams[t][key]) == float:
                    elem.text = "%0.4e" % (tparams[t][key])
                if type(tparams[t][key]) == str:
                    elem.text = "%s" % (tparams[t][key])
    
                    
    mf1comments = SubElement(leapr, 'mf1comments')
    mf1comments.text = comments
    
    # print prettify(model)
    
    if (xml_output == ""):
        return prettify(model)
    else:
        f = open(xml_output, "w")
        f.write(prettify(model))
        f.close()

if __name__ == '__main__':
  if (len(argv)<=2):
    print("Parses a LEAPR input file and creates a XML file")
    print("Usage: %s leapr_input output.xml" % (argv[0]))
    raise SystemExit
  leapr_input = argv[1]
  xml_output = argv[2]
  leapr2xml(leapr_input, xml_output)
 
