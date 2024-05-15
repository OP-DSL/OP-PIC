import os
from math import ceil
from typing import Optional
from util import ABDC, findIdx

from jinja2 import Environment, FileSystemLoader, pass_context

import op as OP

env = Environment(
    loader=FileSystemLoader(os.path.join(os.path.dirname(__file__), "../resources/templates")),
    lstrip_blocks=True,
    trim_blocks=True
)

def direct(x, loop: Optional[OP.Loop] = None) -> bool:
    
    if isinstance(x, (OP.ArgDat, OP.ArgIdx)):
        if x.map_id is None and x.p2c_id is None:
            return True
        else:
            return False

    if isinstance(x, OP.Dat):
        assert loop is not None

        arg = loop.args[x.arg_id]
        assert isinstance(arg, OP.ArgDat)

        if arg.map_id is None and arg.p2c_id is None:
            return True
        else:
            return False

    if isinstance(x, OP.Loop) and len(x.maps) == 0:
        direct = True
        for arg in x.args:
            if isinstance(arg, (OP.ArgDat, OP.ArgIdx)) and (arg.map_id is not None or arg.p2c_id is not None):
                direct = False
        return direct

    return False

def indirect(x, loop: Optional[OP.Loop] = None) -> bool:
    return not direct(x, loop)

def p2c_mapped(x, loop: Optional[OP.Loop] = None) -> bool:
    if isinstance(x, OP.Loop):
        for arg in x.args:
            if isinstance(arg, (OP.ArgDat, OP.ArgIdx)) and arg.p2c_id is not None:
                return True        
    if isinstance(x, (OP.ArgDat, OP.ArgIdx)) and x.p2c_id is not None:
        return True
    if isinstance(x, OP.Map):
        assert loop is not None
        arg = loop.args[x.arg_id]
        assert arg is not None
        return arg.p2c_id is not None
    return False    

def double_indirect(x, loop: Optional[OP.Loop] = None) -> bool:
    if isinstance(x, OP.Loop):
        for arg in x.args:
            if isinstance(arg, (OP.ArgDat, OP.ArgIdx)) and arg.p2c_id is not None and arg.map_id is not None:
                return True        
    if isinstance(x, (OP.ArgDat, OP.ArgIdx)) and x.p2c_id is not None and x.map_id is not None:
        return True   
    return False  

def double_indirect_reduc(x, loop: Optional[OP.Loop] = None) -> bool:
    if x is not None and isinstance(x, OP.Loop):
        for arg in x.args:
            if isinstance(arg, OP.ArgDat) and double_indirect(arg) and arg.access_type in [OP.AccessType.INC, OP.AccessType.MAX, OP.AccessType.MIN]:
                return True
    return False  

def indirect_reduction(x, loop: OP.Loop, check_all: Optional[bool] = False) -> bool:
    if direct(x, loop):
        return False
    
    if not check_all:
        arg = loop.args[x.arg_id]
        assert isinstance(arg, OP.ArgDat)
        if arg.access_type in [OP.AccessType.INC,OP.AccessType.MIN,OP.AccessType.MAX]:
            return True
    return False

def sr_both_direct_and_indirect_reduction(x, loop: OP.Loop) -> bool:
    assert isinstance(x, OP.Dat)
    assert isinstance(loop, OP.Loop)

    direct_mapped = False
    indir_reduc_mapped = False

    for arg in loop.args:
        if isinstance(arg, OP.ArgDat) and arg.dat_id == x.id:
            if arg.access_type in [OP.AccessType.INC,OP.AccessType.MIN,OP.AccessType.MAX]:
                indir_reduc_mapped = True
            else:
                 direct_mapped = True
    return indir_reduc_mapped and direct_mapped

env.tests["injected_loop"] = lambda loop: loop.iterator_type == OP.IterateType.injected
env.tests["particle_loop"] = lambda loop: loop.iterator_set.cell_set != None
env.tests["direct"] = direct
env.tests["indirect"] = indirect
env.tests["p2c_mapped"] = p2c_mapped
env.tests["double_indirect"] = double_indirect
env.tests["double_indirect_reduc"] = double_indirect_reduc
env.tests["indirect_reduction"] = indirect_reduction
env.tests["sr_both_direct_and_indirect_reduction"] = sr_both_direct_and_indirect_reduction

env.tests["map_not_dat"] = lambda dat: dat.flag == False

def soa(x, loop: Optional[OP.Loop] = None) -> bool:
    if isinstance(x, OP.ArgDat):
        assert loop is not None

    if isinstance(x, OP.ArgDat) and loop is not None:
        return loop.dats[x.dat_id].soa

    if isinstance(x, OP.Dat):
        return x.soa

    return False


env.tests["soa"] = soa

env.tests["opt"] = lambda arg, loop=None: hasattr(arg, "opt") and arg.opt

env.tests["dat"] = lambda arg, loop=None: isinstance(arg, OP.ArgDat)
env.tests["gbl"] = lambda arg, loop=None: isinstance(arg, OP.ArgGbl)
env.tests["idx"] = lambda arg, loop=None: isinstance(arg, OP.ArgIdx)
env.tests["info"] = lambda arg, loop=None: isinstance(arg, OP.ArgInfo)

env.tests["vec"] = lambda arg, loop=None: isinstance(arg, OP.ArgDat) and arg.map_idx is not None and arg.map_idx < -1
env.tests["runtime_map_idx"] = (
    lambda arg, loop=None: isinstance(arg, OP.ArgDat) and arg.map_id is not None and arg.map_idx is None
)

env.tests["read"] = lambda arg, loop=None: hasattr(arg, "access_type") and arg.access_type == OP.AccessType.READ
env.tests["write"] = lambda arg, loop=None: hasattr(arg, "access_type") and arg.access_type == OP.AccessType.WRITE
env.tests["read_write"] = lambda arg, loop=None: hasattr(arg, "access_type") and arg.access_type == OP.AccessType.RW

env.tests["inc"] = lambda arg, loop=None: hasattr(arg, "access_type") and arg.access_type == OP.AccessType.INC
env.tests["min"] = lambda arg, loop=None: hasattr(arg, "access_type") and arg.access_type == OP.AccessType.MIN
env.tests["max"] = lambda arg, loop=None: hasattr(arg, "access_type") and arg.access_type == OP.AccessType.MAX

env.tests["work"] = lambda arg, loop=None: hasattr(arg, "access_type") and arg.access_type == OP.AccessType.WORK

env.tests["read_or_write"] = lambda arg, loop=None: hasattr(arg, "access_type") and arg.access_type in [
    OP.AccessType.READ,
    OP.AccessType.WRITE,
    OP.AccessType.RW,
]

env.tests["min_or_max"] = lambda arg, loop=None: hasattr(arg, "access_type") and arg.access_type in [
    OP.AccessType.MIN,
    OP.AccessType.MAX,
]

env.tests["reduction"] = lambda arg, loop=None: hasattr(arg, "access_type") and arg.access_type in [
    OP.AccessType.INC,
    OP.AccessType.MIN,
    OP.AccessType.MAX,
]

def not_already_mapped(arg, loop: Optional[OP.Loop] = None) -> bool:
    if isinstance(arg, OP.ArgDat):
        assert loop is not None
        dat = loop.dats[arg.dat_id]
        for a in loop.args:
            if isinstance(a, OP.ArgDat) and a.dat_id == dat.id and a.id < arg.id and a.access_type in [OP.AccessType.INC,OP.AccessType.MIN,OP.AccessType.MAX]:
                return False
        return True       
    return False   

def same_iter_set_dat(x, loop: OP.Loop) -> bool:
    if isinstance(x, OP.Dat):
        return x.set == loop.iterator_set.name
    elif isinstance(x, OP.ArgDat):
        dat = loop.dats[x.dat_id]
        assert dat is not None
        return dat.set == loop.iterator_set.name
    return False

def not_same_iter_set_dat(x, loop: OP.Loop) -> bool:
    return not same_iter_set_dat(x, loop)

env.tests["not_already_mapped"] = not_already_mapped
env.tests["same_iter_set_dat"] = same_iter_set_dat
env.tests["not_same_iter_set_dat"] = not_same_iter_set_dat

def read_in(dat: OP.Dat, loop: OP.Loop) -> bool:
    for arg in loop.args:
        if not isinstance(arg, OP.ArgDat):
            continue

        if arg.dat_id == dat.id and arg.access_type != OP.AccessType.READ:
            return False

    return True


env.tests["read_in"] = read_in
env.tests["instance"] = lambda x, c: isinstance(x, c)


@pass_context
def select2_filter(context, xs, *tests):
    ys = []
    for x in xs:
        for test in tests:
            if context.environment.call_test(test, x):
                ys.append(x)
                break

    return ys

env.filters["select2"] = select2_filter


def unpack(tup):
    if not isinstance(tup, tuple):
        return tup

    return tup[0]


def test_to_filter(filter_, key=unpack):
    return lambda xs, loop=None: list(filter(lambda x: env.tests[filter_](key(x), loop), xs))  # type: ignore


env.filters["direct"] = test_to_filter("direct")
env.filters["indirect"] = test_to_filter("indirect")
env.filters["not_already_mapped"] = test_to_filter("not_already_mapped")
env.filters["same_iter_set_dat"] = test_to_filter("same_iter_set_dat")
env.filters["not_same_iter_set_dat"] = test_to_filter("not_same_iter_set_dat")

env.filters["soa"] = test_to_filter("soa")

env.filters["opt"] = test_to_filter("opt")

env.filters["dat"] = test_to_filter("dat")
env.filters["gbl"] = test_to_filter("gbl")
env.filters["idx"] = test_to_filter("idx")
env.filters["info"] = test_to_filter("info")

env.filters["vec"] = test_to_filter("vec")
env.filters["runtime_map_idx"] = test_to_filter("runtime_map_idx")

env.filters["read"] = test_to_filter("read")
env.filters["write"] = test_to_filter("write")
env.filters["read_write"] = test_to_filter("read_write")

env.filters["inc"] = test_to_filter("inc")
env.filters["min"] = test_to_filter("min")
env.filters["max"] = test_to_filter("max")

env.filters["read_or_write"] = test_to_filter("read_or_write")

env.filters["reduction"] = test_to_filter("reduction")
env.filters["min_or_max"] = test_to_filter("min_or_max")

env.filters["index"] = lambda xs, x: xs.index(x) if x is not None and x in xs else -1

env.filters["round_up"] = lambda x, b: b * ceil(x / b)


