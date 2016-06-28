#!/usr/bin/python
#testing

import sys
import system_def

para = system_def.SystemDef(sys.argv[1])

print para.pts, para.delta, para.size, para.input_external_field, para.beta


