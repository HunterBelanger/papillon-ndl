#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import sys
import leapr2xml
import leapr_interpolator

if __name__ == '__main__':
  if (len(sys.argv)<=1):
    print("Checks a LEAPR input file")
    print("Usage: %s leapr_input" % (sys.argv[0]))
    raise SystemExit
  leapr_input = sys.argv[1]
  xml_temp = leapr2xml.leapr2xml(leapr_input)
  temp_model = leapr_interpolator.xml2obj(xml_temp)
  leapr_interpolator.check_input(temp_model)

