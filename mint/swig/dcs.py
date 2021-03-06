# This file was automatically generated by SWIG (http://www.swig.org).
# Version 1.3.40
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.
# This file is compatible with both classic and new-style classes.

from sys import version_info
if version_info >= (2,6,0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        fp = None
        try:
            fp, pathname, description = imp.find_module('_dcs', [dirname(__file__)])
        except ImportError:
            import _dcs
            return _dcs
        if fp is not None:
            try:
                _mod = imp.load_module('_dcs', fp, pathname, description)
            finally:
                fp.close()
            return _mod
    _dcs = swig_import_helper()
    del swig_import_helper
else:
    import _dcs
del version_info
try:
    _swig_property = property
except NameError:
    pass # Python < 2.2 doesn't have 'property'.
def _swig_setattr_nondynamic(self,class_type,name,value,static=1):
    if (name == "thisown"): return self.this.own(value)
    if (name == "this"):
        if type(value).__name__ == 'SwigPyObject':
            self.__dict__[name] = value
            return
    method = class_type.__swig_setmethods__.get(name,None)
    if method: return method(self,value)
    if (not static) or hasattr(self,name):
        self.__dict__[name] = value
    else:
        raise AttributeError("You cannot add attributes to %s" % self)

def _swig_setattr(self,class_type,name,value):
    return _swig_setattr_nondynamic(self,class_type,name,value,0)

def _swig_getattr(self,class_type,name):
    if (name == "thisown"): return self.this.own()
    method = class_type.__swig_getmethods__.get(name,None)
    if method: return method(self)
    raise AttributeError(name)

def _swig_repr(self):
    try: strthis = "proxy of " + self.this.__repr__()
    except: strthis = ""
    return "<%s.%s; %s >" % (self.__class__.__module__, self.__class__.__name__, strthis,)

try:
    _object = object
    _newclass = 1
except AttributeError:
    class _object : pass
    _newclass = 0



def version():
  return _dcs.version()
version = _dcs.version

def info():
  return _dcs.info()
info = _dcs.info
class MParameters(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, MParameters, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, MParameters, name)
    __repr__ = _swig_repr
    def __init__(self, *args): 
        this = _dcs.new_MParameters(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_setmethods__["energy"] = _dcs.MParameters_energy_set
    __swig_getmethods__["energy"] = _dcs.MParameters_energy_get
    if _newclass:energy = _swig_property(_dcs.MParameters_energy_get, _dcs.MParameters_energy_set)
    __swig_destroy__ = _dcs.delete_MParameters
    __del__ = lambda self : None;
MParameters_swigregister = _dcs.MParameters_swigregister
MParameters_swigregister(MParameters)
cvar = _dcs.cvar

class Device(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, Device, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, Device, name)
    __repr__ = _swig_repr
    __swig_setmethods__["id"] = _dcs.Device_id_set
    __swig_getmethods__["id"] = _dcs.Device_id_get
    if _newclass:id = _swig_property(_dcs.Device_id_get, _dcs.Device_id_set)
    __swig_setmethods__["channel_z"] = _dcs.Device_channel_z_set
    __swig_getmethods__["channel_z"] = _dcs.Device_channel_z_get
    if _newclass:channel_z = _swig_property(_dcs.Device_channel_z_get, _dcs.Device_channel_z_set)
    __swig_setmethods__["z_pos"] = _dcs.Device_z_pos_set
    __swig_getmethods__["z_pos"] = _dcs.Device_z_pos_get
    if _newclass:z_pos = _swig_property(_dcs.Device_z_pos_get, _dcs.Device_z_pos_set)
    def __init__(self, *args): 
        this = _dcs.new_Device(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _dcs.delete_Device
    __del__ = lambda self : None;
Device_swigregister = _dcs.Device_swigregister
Device_swigregister(Device)

class BPM(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, BPM, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, BPM, name)
    __repr__ = _swig_repr
    __swig_setmethods__["id"] = _dcs.BPM_id_set
    __swig_getmethods__["id"] = _dcs.BPM_id_get
    if _newclass:id = _swig_property(_dcs.BPM_id_get, _dcs.BPM_id_set)
    __swig_setmethods__["channel_x"] = _dcs.BPM_channel_x_set
    __swig_getmethods__["channel_x"] = _dcs.BPM_channel_x_get
    if _newclass:channel_x = _swig_property(_dcs.BPM_channel_x_get, _dcs.BPM_channel_x_set)
    __swig_setmethods__["channel_y"] = _dcs.BPM_channel_y_set
    __swig_getmethods__["channel_y"] = _dcs.BPM_channel_y_get
    if _newclass:channel_y = _swig_property(_dcs.BPM_channel_y_get, _dcs.BPM_channel_y_set)
    __swig_setmethods__["channel_z"] = _dcs.BPM_channel_z_set
    __swig_getmethods__["channel_z"] = _dcs.BPM_channel_z_get
    if _newclass:channel_z = _swig_property(_dcs.BPM_channel_z_get, _dcs.BPM_channel_z_set)
    __swig_setmethods__["x"] = _dcs.BPM_x_set
    __swig_getmethods__["x"] = _dcs.BPM_x_get
    if _newclass:x = _swig_property(_dcs.BPM_x_get, _dcs.BPM_x_set)
    __swig_setmethods__["y"] = _dcs.BPM_y_set
    __swig_getmethods__["y"] = _dcs.BPM_y_get
    if _newclass:y = _swig_property(_dcs.BPM_y_get, _dcs.BPM_y_set)
    __swig_setmethods__["z_pos"] = _dcs.BPM_z_pos_set
    __swig_getmethods__["z_pos"] = _dcs.BPM_z_pos_get
    if _newclass:z_pos = _swig_property(_dcs.BPM_z_pos_get, _dcs.BPM_z_pos_set)
    def __init__(self, *args): 
        this = _dcs.new_BPM(*args)
        try: self.this.append(this)
        except: self.this = this
    def info(self): return _dcs.BPM_info(self)
    __swig_destroy__ = _dcs.delete_BPM
    __del__ = lambda self : None;
BPM_swigregister = _dcs.BPM_swigregister
BPM_swigregister(BPM)

class Orbit(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, Orbit, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, Orbit, name)
    __repr__ = _swig_repr
    def __init__(self): 
        this = _dcs.new_Orbit()
        try: self.this.append(this)
        except: self.this = this
    __swig_setmethods__["bpms"] = _dcs.Orbit_bpms_set
    __swig_getmethods__["bpms"] = _dcs.Orbit_bpms_get
    if _newclass:bpms = _swig_property(_dcs.Orbit_bpms_get, _dcs.Orbit_bpms_set)
    def get(self, *args): return _dcs.Orbit_get(self, *args)
    def __getitem__(self, *args): return _dcs.Orbit___getitem__(self, *args)
    __swig_destroy__ = _dcs.delete_Orbit
    __del__ = lambda self : None;
Orbit_swigregister = _dcs.Orbit_swigregister
Orbit_swigregister(Orbit)

class Func_1d(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, Func_1d, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, Func_1d, name)
    __repr__ = _swig_repr
    def __init__(self, *args): 
        this = _dcs.new_Func_1d(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_setmethods__["f"] = _dcs.Func_1d_f_set
    __swig_getmethods__["f"] = _dcs.Func_1d_f_get
    if _newclass:f = _swig_property(_dcs.Func_1d_f_get, _dcs.Func_1d_f_set)
    __swig_setmethods__["n"] = _dcs.Func_1d_n_set
    __swig_getmethods__["n"] = _dcs.Func_1d_n_get
    if _newclass:n = _swig_property(_dcs.Func_1d_n_get, _dcs.Func_1d_n_set)
    def get(self, *args): return _dcs.Func_1d_get(self, *args)
    def set(self, *args): return _dcs.Func_1d_set(self, *args)
    def sum(self): return _dcs.Func_1d_sum(self)
    def __getitem__(self, *args): return _dcs.Func_1d___getitem__(self, *args)
    def __setitem__(self, *args): return _dcs.Func_1d___setitem__(self, *args)
    def __len__(self): return _dcs.Func_1d___len__(self)
    def __str__(self): return _dcs.Func_1d___str__(self)
    def __add__(self, *args): return _dcs.Func_1d___add__(self, *args)
    __swig_destroy__ = _dcs.delete_Func_1d
    __del__ = lambda self : None;
Func_1d_swigregister = _dcs.Func_1d_swigregister
Func_1d_swigregister(Func_1d)


def get_device_info(*args):
  return _dcs.get_device_info(*args)
get_device_info = _dcs.get_device_info

def get_device_val(*args):
  return _dcs.get_device_val(*args)
get_device_val = _dcs.get_device_val

def set_device_val(*args):
  return _dcs.set_device_val(*args)
set_device_val = _dcs.set_device_val

def get_device_td(*args):
  return _dcs.get_device_td(*args)
get_device_td = _dcs.get_device_td

def test_func_1d(*args):
  return _dcs.test_func_1d(*args)
test_func_1d = _dcs.test_func_1d

def test_func_1d_2(*args):
  return _dcs.test_func_1d_2(*args)
test_func_1d_2 = _dcs.test_func_1d_2

def get_parameters():
  return _dcs.get_parameters()
get_parameters = _dcs.get_parameters

def get_bpm_val(*args):
  return _dcs.get_bpm_val(*args)
get_bpm_val = _dcs.get_bpm_val

def get_orbit():
  return _dcs.get_orbit()
get_orbit = _dcs.get_orbit


