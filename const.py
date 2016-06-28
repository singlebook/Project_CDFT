#!/usr/bin/python
#In Python, any variable can be re-bound at will -- and modules don't let you define special methods such as an 
#instance's __setattr__ to stop attribute re-binding. Easy solution (in Python 2.1 and up): use an instance as  "module"
#
#http://code.activestate.com/recipes/65207-constants-in-python/
#
#Filename: const.py

class _const:
    class ConstError(TypeError):
        pass

    def __setattr__(self, name, value):
        if self.__dict__.has_key(name):
            raise self.ConstError, "Can't rebind const instance attribute (%s)" % name

        self.__dict__[name] = value
        
    def __delattr__(self, name):
		if self.__dict__.has_key(name):
			raise self.ConstError, "Can't unbind const const instance attribute (%s)" % name
		
		raise AttributeError, "const instance has no attribute '%s'" % name


import sys
sys.modules[__name__] = _const()

