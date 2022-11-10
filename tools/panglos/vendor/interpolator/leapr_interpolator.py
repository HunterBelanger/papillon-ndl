#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import re
import xml.sax.handler
import numpy as np
import sys
import copy
import leapr2xml

def xml2obj(src):
    """
    A simple function to converts XML data into native Python object.
    """

    non_id_char = re.compile('[^_0-9a-zA-Z]')
    def _name_mangle(name):
        return non_id_char.sub('_', name)

    class DataNode(object):
        def __init__(self):
            self._attrs = {}    # XML attributes and child elements
            self.data = None    # child text data
        def __len__(self):
            # treat single element as a list of 1
            return 1
        def __getitem__(self, key):
            if isinstance(key, str):
                return self._attrs.get(key,None)
            else:
                return [self][key]
        def __contains__(self, name):
            return name in self._attrs
        def __nonzero__(self):
            return bool(self._attrs or self.data)
        def __getattr__(self, name):
            if name.startswith('__'):
                # need to do this for Python special methods???
                raise AttributeError(name)
            return self._attrs.get(name,None)
        def _add_xml_attr(self, name, value):
            if name in self._attrs:
                # multiple attribute of the same name are represented by a list
                children = self._attrs[name]
                if not isinstance(children, list):
                    children = [children]
                    self._attrs[name] = children
                children.append(value)
            else:
                self._attrs[name] = value
        def __str__(self):
            return self.data or ''
        def __repr__(self):
            items = sorted(self._attrs.items())
            if self.data:
                items.append(('data', self.data))
            return u'{%s}' % ', '.join([u'%s:%s' % (k,repr(v)) for k,v in items])

    class TreeBuilder(xml.sax.handler.ContentHandler):
        def __init__(self):
            self.stack = []
            self.root = DataNode()
            self.current = self.root
            self.text_parts = []
        def startElement(self, name, attrs):
            self.stack.append((self.current, self.text_parts))
            self.current = DataNode()
            self.text_parts = []
            # xml attributes --> python attributes
            for k, v in attrs.items():
                self.current._add_xml_attr(_name_mangle(k), v)
        def endElement(self, name):
            text = ''.join(self.text_parts).strip()
            if text:
                self.current.data = text
            if self.current._attrs:
                obj = self.current
            else:
                # a text only node is simply represented by the string
                obj = text or ''
            self.current, self.text_parts = self.stack.pop()
            self.current._add_xml_attr(_name_mangle(name), obj)
        def characters(self, content):
            self.text_parts.append(content)

    builder = TreeBuilder()
    if isinstance(src,str):
        xml.sax.parseString(src, builder)
    else:
        xml.sax.parse(src, builder)
    return list(builder.root._attrs.values())[0]

def interpolate(x, x1, y1, x2, y2, mode="linear"):
    """
    Interpolates between two points in different ways
    """
    if mode == "linear":
        return (y2-y1)/(x2-x1)*(x-x1)+y1
    if mode == "arrhenius":
        return np.exp((np.log(y2)-np.log(y1))/(1.0/x2-1.0/x1)*(1.0/x-1.0/x1)+np.log(y1))
    if mode == "reciprocal":
        return 1.0/((1.0/y2-1.0/y1)/(x2-x1)*(x-x1)+1.0/y1)
        
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
    
def new_temperature(model, new_t):
    """
    Interpolates a model tree for a new temperature
    """
    temps = [float(t.temp) for t in model.leapr.temperature]

    if (model.leapr.nsk == None):
        nsk = 0
    else:
        nsk = int(model.leapr.nsk)
    if ((new_t >= min(temps)) and (new_t <= max(temps))):
        for t in model.leapr.temperature:
            if (new_t == float(t.temp)):
                print("Temperature already in the input file")
                return t
        deltat = np.Inf
        for t in model.leapr.temperature:
            if (((new_t - float(t.temp)) > 0.0) and (abs(new_t - float(t.temp)) < deltat)):
                t1 = t
                deltat = new_t - float(t.temp)
        deltat = np.Inf
        for t in model.leapr.temperature:
            if (((new_t - float(t.temp)) < 0.0) and (abs(new_t - float(t.temp)) < deltat)):
                t2 = t
                deltat = new_t - float(t.temp)
        print("Interpolating %0.2f between %0.2f and %0.2f" % (new_t, float(t1.temp), float(t2.temp)))
    elif (new_t < min(temps)):
        new_t1 = np.Inf
        for t in model.leapr.temperature:
            if (float(t.temp) < new_t1):
                t1 = t
                new_t1 = float(t.temp)
        new_t2 = np.Inf
        for t in model.leapr.temperature:
            if ((new_t1 < float(t.temp)) and (float(t.temp) < new_t2)):
                t2 = t
                new_t2 = float(t.temp)
        print("Extrapolating %0.2f from %0.2f and %0.2f" % (new_t, float(t1.temp), float(t2.temp)))
    elif (new_t > max(temps)):
        new_t1 = 0.0
        for t in model.leapr.temperature:
            if (new_t1 < float(t.temp)):
                t1 = t
                new_t1 = float(t.temp)
        new_t2 = 0.0
        for t in model.leapr.temperature:
            if ((new_t1 > float(t.temp)) and (float(t.temp) > new_t2)):
                t2 = t
                new_t2 = float(t.temp)
        print("Extrapolating %0.2f from %0.2f and %0.2f" % (new_t, float(t1.temp), float(t2.temp)))
    temp1 = float(t1.temp)
    temp2 = float(t2.temp)
    if (temp1 == temp2):
        print("Need at least two different temperatures for interpolation")
        return None
#
# Continuous frequency spectrum
#
    delta1 = float(t1.delta)
    delta2 = float(t2.delta)
    new_delta = interpolate(new_t, temp1, delta1, temp2, delta2, mode="linear") 
    ni1 = float(t1.ni)
    ni2 = float(t2.ni)
    rho1 = np.fromstring(t1.rho, sep= " ")
    rho2 = np.fromstring(t2.rho, sep= " ")
    erg1 = np.linspace(0.0, ni1*delta1, int(ni1))
    erg2 = np.linspace(0.0, ni2*delta2, int(ni2))
    new_erg = np.arange(0.0, min(max(erg1),max(erg2)),new_delta)
    new_ni = len(new_erg)
    rho1_int=np.interp(new_erg, erg1, rho1)
    rho2_int=np.interp(new_erg, erg2, rho2)
    new_rho = interpolate(new_t, temp1, rho1_int, temp2, rho2_int, mode="linear")
#
# Translational properties and continuous frequency spectrum normalization
#
    twt1 = float(t1.twt)
    twt2 = float(t2.twt)
    new_twt = interpolate(new_t, temp1, twt1, temp2, twt2, mode="linear") 
    c1 = float(t1.c)
    c2 = float(t2.c)
    if (c1*c2 == 0.0):
        new_c = 0.0
    else:
        new_c = 1.0/new_twt*interpolate(new_t, temp1, c1*twt1, temp2, c2*twt2, mode="arrhenius") 
    tbeta1 = float(t1.tbeta)
    tbeta2 = float(t2.tbeta)
    new_tbeta = interpolate(new_t, temp1, tbeta1, temp2, tbeta2, mode="linear") 
#
# Oscillators
#
    if (int(t1.nd) != int(t2.nd)):
        print("Different number of oscillators at different temperatures")
        return None
    else:
        new_nd = int(t1.nd)
    bdel1 = np.fromstring(t1.bdel, sep= " ")
    bdel2 = np.fromstring(t2.bdel, sep= " ")
    new_bdel = interpolate(new_t, temp1, bdel1, temp2, bdel2, mode="linear") 
    adel1 = np.fromstring(t1.adel, sep= " ")
    adel2 = np.fromstring(t2.adel, sep= " ")
    new_adel = interpolate(new_t, temp1, adel1, temp2, adel2, mode="linear") 
#
# Structure correction factors
#
    if (nsk > 0):
        dka1 = float(t1.dka)
        dka2 = float(t2.dka)
        new_dka = interpolate(new_t, temp1, dka1, temp2, dka2, mode="linear") 
        nka1 = float(t1.nka)
        nka2 = float(t2.nka)
        skappa1 = np.fromstring(t1.skappa, sep= " ")
        skappa2 = np.fromstring(t2.skappa, sep= " ")
        kappa1 = np.linspace(0.0, nka1*dka1,nka1)
        kappa2 = np.linspace(0.0, nka2*dka2,nka2)
        new_kappa = np.arange(0.0, min(max(kappa1),max(kappa2)),new_dka)
        new_nka = len(new_kappa)
        skappa1_int=np.interp(new_kappa, kappa1, skappa1)
        skappa2_int=np.interp(new_kappa, kappa2, skappa2)
        new_skappa = interpolate(new_t, temp1, skappa1_int, temp2, skappa2_int, mode="linear")
#
# Save data
#
    new_temperature = copy.copy(model.leapr.temperature[0])
    new_temperature.temp = "%0.2f" % (new_t)
    new_temperature.delta = "%0.4e" % (new_delta)
    new_temperature.ni = "%i" % (new_ni)
#    s = np.array2string(new_rho,  max_line_width=55, formatter={'float_kind':lambda x: "%.4e" % x})
#    s = re.sub('[\[\]]', '', s)
#    s = re.sub('\n ', ' \n', s)
    s = five_columns(new_rho)
    new_temperature.rho = s
    new_temperature.twt = "%0.4e" % (new_twt)
    new_temperature.c = "%0.4e" % (new_c)
    new_temperature.tbeta = "%0.4e" % (new_tbeta)
    new_temperature.nd = "%i" % (new_nd)
#    s = np.array2string(new_bdel,  max_line_width=55, formatter={'float_kind':lambda x: "%.4e" % x})
#    s = re.sub('[\[\]]', '', s)
#    s = re.sub('\n ', ' \n', s)
    s = five_columns(new_bdel)
    new_temperature.bdel = s
#    s = np.array2string(new_adel,  max_line_width=55, formatter={'float_kind':lambda x: "%.4e" % x})
#    s = re.sub('[\[\]]', '', s)
#    s = re.sub('\n ', ' \n', s)
    s = five_columns(new_adel)
    new_temperature.adel = s
    if (nsk > 0):
        new_temperature.dka = "%0.4e" % (new_dka)
        new_temperature.nka = "%i" % (new_nka)
#        s = np.array2string(new_skappa,  max_line_width=55, formatter={'float_kind':lambda x: "%.4e" % x})
#        s = re.sub('[\[\]]', '', s)
#        s = re.sub('\n ', ' \n', s)
        s = five_columns(new_skappa)
        new_temperature.skappa = s
    return new_temperature

def print_input(model):
    print(generate_input(model))

def generate_input(model):
    """
    Receives a model tree and returns a LEAPR input as a string
    """

    inp = ""
    inp = inp + "leapr\n"
    inp = inp + model.leapr.nout + " / " + model.leapr.c1comment + "\n"
    inp = inp + model.leapr.title + " / " + model.leapr.c2comment + "\n"
    inp = inp + model.leapr.ntempr + " " + model.leapr.iprint + " " + model.leapr.nphon + " / " + " " + model.leapr.c3comment + "\n"
    inp = inp + model.leapr.mat + " " + model.leapr.za + " " + model.leapr.isabt + " " + model.leapr.ilog + " / " + model.leapr.c4comment + "\n"
    inp = inp + model.leapr.awr + " " + model.leapr.spr + " " + model.leapr.npr + " " + model.leapr.iel + " " + model.leapr.ncold + " " + model.leapr.nsk + " / " + model.leapr.c5comment + "\n"
    inp = inp + model.leapr.nss + " " + model.leapr.b7 + " " + model.leapr.aws + " " + model.leapr.sps + " " + model.leapr.mss + " / " + model.leapr.c6comment + "\n"
    inp = inp + model.leapr.nalpha + " " + model.leapr.nbeta + " " + model.leapr.lat + " / " + model.leapr.c7comment + "\n"
    inp = inp + model.leapr.alphas + " / " + model.leapr.c8comment + "\n"
    inp = inp + model.leapr.betas + " / " + model.leapr.c9comment + "\n"
    for t in model.leapr.temperature:
        inp = inp + t.temp + " / " + t.c10comment + "\n"
        inp = inp + t.delta + " " + t.ni + " / " + t.c11comment + "\n"
        inp = inp + t.rho + " / " + t.c12comment + "\n"
        inp = inp + t.twt + " " + t.c + " " + t.tbeta + " / " + t.c13comment + "\n"
        inp = inp + t.nd + " / " + t.c14comment + "\n"
        if (int(t.nd) > 0):
            inp = inp + t.bdel + " / " + t.c15comment + "\n"
            inp = inp + t.adel + " / " + t.c16comment + "\n"
        if (int(model.leapr.nsk) > 0):
            inp = inp + t.nka + " " + t.dka + " / " + t.c17comment + "\n"
            inp = inp + t.skappa + " / " + t.c18comment + "\n"
            inp = inp + t.cfrac + " / " + t.c19comment + "\n"
    inp = inp + model.leapr.mf1comments + "\n"
    inp = inp + "/\n"
    inp = inp + "stop\n"
    return inp

def check_input(model):
    """
    Checks a model tree
    """

    eps = 1e-4
    print("Checking LEAPR input")
    error = 0
    warning = 0
    for tag in [model.leapr.c1comment, model.leapr.c2comment, model.leapr.c3comment, model.leapr.c4comment, model.leapr.c5comment, \
                model.leapr.c6comment, model.leapr.c7comment, model.leapr.c8comment, model.leapr.c9comment, \
                model.leapr.nout, model.leapr.title, model.leapr.ntempr, model.leapr.iprint, model.leapr.nphon, \
                model.leapr.mat, model.leapr.za, model.leapr.isabt, model.leapr.ilog, \
                model.leapr.awr, model.leapr.spr, model.leapr.npr, model.leapr.iel, model.leapr.ncold, model.leapr.nsk, \
                model.leapr.nss, model.leapr.b7, model.leapr.aws, model.leapr.sps, model.leapr.mss, \
                model.leapr.nalpha, model.leapr.nbeta, model.leapr.lat, \
                model.leapr.alphas, model.leapr.betas, \
                model.leapr.mf1comments]:
        if type(tag)==list:
            print("Error: some of the tags appear to be repeated")
            error = error + 1
#
# Card 1
#
    if ((model.leapr.nout == None) or not (int(model.leapr.nout) > 0)):
        print("Error: nout required")
        error = error + 1 
#
# Card 2
#
    if (model.leapr.title == None):
        model.leapr.title = ""
    print("TITLE: %s" % (model.leapr.title))
#
# Card 3
#
    if ((model.leapr.ntempr == None) or not (int(model.leapr.ntempr) > 0)):
        print("Error: ntempr required")
        error = error + 1 
    if (int(model.leapr.ntempr) != len(model.leapr.temperature)):
        print("Error: ntempr not equal to number of temperatures")
        error = error + 1 
    if (model.leapr.iprint == None):
        model.leapr.iprint = ""
        if (model.leapr.nphon != None):
            print("Warning: using default 0 value for iprint")
            warning = warning + 1
            model.leapr.iprint = "0"
    if (model.leapr.nphon == None):
        model.leapr.nphon = ""
#
# Card 4
#
    if (model.leapr.mat == None):
        print("Error: mat required")
        error = error + 1 
    if (model.leapr.za == None):
        print("Error: za required")
        error = error + 1 
    if (model.leapr.isabt == None):
        model.leapr.isabt = ""
        if (model.leapr.ilog != None):
            print("Warning: using default 0 value for isabt")
            warning = warning + 1
            model.leapr.isabt = "0"
    if (model.leapr.ilog == None):
        model.leapr.ilog = ""
#
# Card 5
#
    if (model.leapr.awr == None):
        print("Error: awr required")
        error = error + 1 
    if (model.leapr.spr == None):
        print("Error: spr required")
        error = error + 1 
    if (model.leapr.npr == None):
        print("Error: npr required")
        error = error + 1 
    if (model.leapr.iel == None):
        model.leapr.iel = ""
        if ((model.leapr.ncold != None) or (model.leapr.nsk != None)):
            print("Warning: using default 0 value for iel")
            warning = warning + 1
            model.leapr.iel = "0"
    if (model.leapr.ncold == None):
        model.leapr.ncold = ""
        if (model.leapr.nsk != None):
            print("Warning: using default 0 value for ncold")
            warning = warning + 1
            model.leapr.ncold = "0"
    if (model.leapr.nsk == None):
        model.leapr.nsk = ""
#
# Card 6
#
    if (model.leapr.nss == None):
        print("Warning: using default 0 value for nss")
        warning = warning + 1
        model.leapr.nss = "0"
    if (int(model.leapr.nss) !=0 and ((model.leapr.b7 == None) or (model.leapr.aws == None) or (model.leapr.sps == None) or (model.leapr.mss == None))):
        print("Error: b7, aws, sps and mss are needed for secondary scatterers")
        error = error + 1 
#
# Card 7
#
    if ((model.leapr.nalpha == None) or not (int(model.leapr.nalpha) > 0)):
        print("Error: nalpha required")
        error = error + 1 
    if ((model.leapr.nbeta == None) or not (int(model.leapr.nbeta) > 0)):
        print("Error: nbeta required")
        error = error + 1 
    if (model.leapr.lat == None):
        model.leapr.lat == ""
    if (int(model.leapr.lat) not in [0, 1]):
        print("Error: lat must be 0 or 1")
        error = error + 1
#
# Card 8
#
    if (model.leapr.alphas == None):
        print("Error: alpha grid required")
        error = error + 1 
    alphas = np.fromstring(model.leapr.alphas, sep= " ")
    if (len(alphas) != int(model.leapr.nalpha)):
        print("Error: alpha grid not nalpha long")
        error = error + 1 
    if (alphas <= 0.0).any():
        print("Error: all alphas must be greater than zero")
        error = error + 1 
    if (not np.all(np.diff(alphas) > 0)):
        print("Error: alphas must be ordered")
        error = error + 1 

#
# Card 9
#
    if (model.leapr.betas == None):
        print("Error: beta grid required")
        error = error + 1 
    betas = np.fromstring(model.leapr.betas, sep= " ")
    if (len(betas) != int(model.leapr.nbeta)):
        print("Error: beta grid not nbeta long")
        error = error + 1 
    if (betas < 0.0).any():
        print("Error: all betas must be greater or equal than zero")
        error = error + 1 
    if (not np.all(np.diff(betas) > 0)):
        print("Error: betas must be ordered")
        error = error + 1 
#
# Temperature loop
#
    for t in model.leapr.temperature:
#
# Card 10
#
        for tag in [t.c10comment, t.c11comment, t.c12comment, t.c13comment, t.c14comment, \
                    t.c15comment, t.c16comment, t.c17comment, t.c18comment, t.c19comment, \
                    t.temp, t.delta, t.ni, t.rho, t.twt, t.c, t.beta, t.nd, t.adel, t.bdel, \
                    t.nka, t.dka, t.skappa, t.cfrac]:
            if type(tag)==list:
                print("Error: some of the tags appear to be repeated")
                error = error + 1
        if (t.temp == None):
            print("Error: temperature required")
            error = error + 1 

        if (float(t.temp) > 0.0):
            if (t.delta == None):
                print("Error: delta required")
                error = error + 1 
            if (t.ni == None):
                print("Error: ni required")
                error = error + 1 
            if (t.rho == None):
                print("Error: rho required")
                error = error + 1 
            rho = np.fromstring(t.rho, sep= " ")
            if (len(rho) != int(t.ni)):
                print("Error: rho not ni long")
                error = error + 1 
            if ((rho < 0.0).any()):
                print("Error: some of the elements of rho are less than zero")
                error = error + 1 
            if (t.twt == None):
                print("Error: twt required")
                error = error + 1 
            if (t.c == None):
                print("Error: c required")
                error = error + 1 
            if (t.tbeta == None):
                print("Error: tbeta required")
                error = error + 1 
            if ((float(t.twt) > 0.0) and (int(model.leapr.iel)>0)):
                print("Error: translational model in an elastic material")
                error = error + 1 
            if (float(t.twt) < 0.0):
                print("Error: twt less than 0.0")
                error = error + 1 
            if (float(t.c) < 0.0):
                print("Error: c less than 0.0")
                error = error + 1 
            if (float(t.tbeta) < 0.0):
                print("Error: tbeta less than 0.0")
                error = error + 1 
            if (t.nd == None):
                print("Error: nd required, set to 0")
                error = error + 1 
                t.nd = "0"
            if (int(t.nd)>0):
                if (t.bdel == None):
                    print("Error: bdel required if nd > 0")
                    error = error + 1 
                bdel = np.fromstring(t.bdel, sep= " ")
                if (len(bdel) != int(t.nd)):
                    print("Error: bdel not nd long")
                    error = error + 1 
                if (t.adel == None):
                    print("Error: adel required if nd > 0")
                    error = error + 1 
                adel = np.fromstring(t.adel, sep= " ")
                tot_weight = float(t.twt)+float(t.tbeta) + np.sum(np.array(adel))
                if (abs(tot_weight - 1.0) > eps):
                    print("Error: sum of weights: %0.4f for temperature %0.2f" % (tot_weight, float(t.temp)))
                    error = error + 1
                if (len(adel) != int(t.nd)):
                    print("Error: adel not nd long")
                    error = error + 1 
            else:
                tot_weight = float(t.twt)+float(t.tbeta)
                if (abs(tot_weight - 1.0) > eps):
                    print("Error: sum of weights: %0.4f for temperature %0.2f" % (tot_weight, float(t.temp)))
                    error = error + 1
                if (t.bdel != None):
                    print("Error: bdel ignored because nd == 0")
                    error = error + 1 
                if (t.adel != None):
                    print("Error: adel ignored because nd == 0")
                    error = error + 1 
            if (int(model.leapr.nsk)>0):
                if (t.nka == None):
                    print("Error: nka required if nsk > 0")
                    error = error + 1 
                if not (int(t.nka) >0):
                    print("Error: nka must be greater than 0")
                    error = error + 1 
                if (t.dka == None):
                    print("Error: dka required if nsk > 0")
                    error = error + 1 
                if not (float(t.dka) >0.0):
                    print("Error: dka must be greater than 0.0")
                    error = error + 1 
                if (t.skappa == None):
                    print("Error: skappa required if nsk > 0")
                    error = error + 1 
                skappa = np.fromstring(t.skappa, sep= " ")
                if (len(skappa) != int(t.nka)):
                    print("Error: skappa not nka long")
                    error = error + 1 
                if ((skappa <= 0.0).any()):
                    print("Error: some of the elements of skappa are less or equal than zero")
                    error = error + 1 
                if (t.cfrac == None):
                    if ((t.ncold == "") or (t.ncold == "0")) :
                        print("Error: cfrac required if nsk > 0 and ncold == 0")
                        error = error + 1 
                elif not (float(t.cfrac) >0.0):
                    if ((t.ncold == "") or (t.ncold == "0")) :
                        print("Error: cfrac must be greater than 0.0 if nsk > 0 and ncold == 0")
                        error = error + 1 
#
# MF=1 Comments
#
    if (model.leapr.mf1comments == None):
        model.leapr.mf1comments = ""
#
# Process comments, if any
#
    if (model.leapr.c1comment == None):
        model.leapr.c1comment = ""
    if (model.leapr.c2comment == None):
        model.leapr.c2comment = ""
    if (model.leapr.c3comment == None):
        model.leapr.c3comment = ""
    if (model.leapr.c4comment == None):
        model.leapr.c4comment = ""
    if (model.leapr.c5comment == None):
        model.leapr.c5comment = ""
    if (model.leapr.c6comment == None):
        model.leapr.c6comment = ""
    if (model.leapr.c7comment == None):
        model.leapr.c7comment = ""
    if (model.leapr.c8comment == None):
        model.leapr.c8comment = ""
    if (model.leapr.c9comment == None):
        model.leapr.c9comment = ""
    for i in range(len(model.leapr.temperature)):
        if (model.leapr.temperature[i].c10comment == None):
            model.leapr.temperature[i].c10comment = ""
        if (model.leapr.temperature[i].c11comment == None):
            model.leapr.temperature[i].c11comment = ""
        if (model.leapr.temperature[i].c12comment == None):
            model.leapr.temperature[i].c12comment = ""
        if (model.leapr.temperature[i].c13comment == None):
            model.leapr.temperature[i].c13comment = ""
        if (model.leapr.temperature[i].c14comment == None):
            model.leapr.temperature[i].c14comment = ""
        if (model.leapr.temperature[i].c15comment == None):
            model.leapr.temperature[i].c15comment = ""
        if (model.leapr.temperature[i].c16comment == None):
            model.leapr.temperature[i].c16comment = ""
        if (model.leapr.temperature[i].c17comment == None):
            model.leapr.temperature[i].c17comment = ""
        if (model.leapr.temperature[i].c18comment == None):
            model.leapr.temperature[i].c18comment = ""
        if (model.leapr.temperature[i].c19comment == None):
            model.leapr.temperature[i].c19comment = ""
#
# Finish
#
    print("Input processed with", error, "errors and", warning, "warnings.")

def leapr_interpolator(base_xml_model, new_leapr_input, t):
    try:
        f = open(base_xml_model)
        f.close()
    except IOError:
        print("Error: %s does not exist" %(sys.argv[1]))
        raise SystemExit

    xmldoc = open(base_xml_model, "r")
    model = xml2obj(xmldoc)

    check_input(model)

    new_temp = new_temperature(model, t)
    if (new_temp != None):
        ntempr = 1
        model.leapr.ntempr = "%i" % (ntempr)
        model.leapr.title = "Model interpolated at %0.2f" % float(new_temp.temp)
        model.leapr.temperature = [new_temp]

        leapr_input = generate_input(model)

        f = open(new_leapr_input, "w")
        f.write(leapr_input)
        f.close()
    
        xml_temp = leapr2xml.leapr2xml(new_leapr_input)
        temp_model = xml2obj(xml_temp)
        check_input(temp_model)
    else:
        print("Error during interpolation")
        raise SystemExit

def getnout(base_xml_model):
    try:
        f = open(base_xml_model)
        f.close()
    except IOError:
        print("Error: %s does not exist" %(sys.argv[1]))
        raise SystemExit

    xmldoc = open(base_xml_model, "r")
    model = xml2obj(xmldoc)
    
    return model.leapr.nout
    
if __name__ == '__main__':
    if (len(sys.argv)<=3):
        print("Usage: %s base_model.xml new_input.leapr temp" % (sys.argv[0]))
        raise SystemExit
    new_leapr_input = sys.argv[2]
    base_xml_model = sys.argv[1]
    try:
        t = float(sys.argv[3])
    except ValueError:
        print("Error: %s does not seem to be a float" %(sys.argv[3]))
        raise SystemExit
    leapr_interpolator(base_xml_model, new_leapr_input, t)


