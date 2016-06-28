#!/usr/bin/python
#testing

import sys
import system_def
import init_sys

para = system_def.SystemDef(sys.argv[1])

init = init_sys.Init(para.input_external_field, para.beta)

print init.field_external
print init.density
