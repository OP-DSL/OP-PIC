import re

from store import Program
from util import SourceBuffer
import op as OP

def populate_arg(arg, mlh: OP.Loop) -> str:
    id1 = ""
    id2 = ""
    id3 = ""
    if isinstance(arg, OP.ArgDat):
        if arg.p2c_id is not None:
            id3 = f"{arg.p2c_id}, "
        if arg.map_id is not None:
            id2 = f"{mlh.map(arg).ptr}, "
        if arg.map_idx >= 0:
            id1 = f"{arg.map_idx}, "
        arg_str = f"opp_arg_dat({mlh.dat(arg).ptr}, {id1}{id2}{id3}OPP_{arg.access_type.name})"
    elif isinstance(arg, OP.ArgGbl):
        arg_str = f'opp_arg_gbl({arg.ptr}, {arg.dim}, "{arg.typ.__repr__()}", OPP_{arg.access_type.name})'
    return arg_str

def get_dh_init_args(move_loop : OP.Loop):
    if move_loop is None:
        return ""
    dh_args = ", " + move_loop.c2c_map.ptr + ", " + move_loop.p2c_map.ptr + ",\n\t\t\t"
    dh_args += ",\n\t\t\t".join(populate_arg(arg, move_loop) for arg in move_loop.args)
    return dh_args

# Augment source program to use generated kernel hosts
def translateProgram(source: str, program: Program, force_soa: bool) -> str:
    buffer = SourceBuffer(source)

    # 1. Remove the defined constants : OPP_* CONST_*
    pattern = re.compile(r'^\s*(OPP_REAL|OPP_BOOL|OPP_INT)\sCONST_.*$')
    lines_removed = 0
    index = buffer.search(pattern)
    while index is not None:
        buffer.remove(index)
        source = buffer.translate()
        buffer = SourceBuffer(source)
        lines_removed += 1
        index = buffer.search(pattern)

    # 2. Update loop calls
    for loop in program.loops:
        if loop.loop_type == OP.LoopType.PAR_LOOP:
            before, after = buffer.get(loop.loc.line - lines_removed - 1).split("opp_par_loop", 1)
            after = re.sub(
                rf"{loop.kernel}\s*,\s*\"{loop.loop_name}\"\s*,\s*", "", after, count=1
            )  # TODO: This assumes that the kernel arg is on the same line as the call
            buffer.update(loop.loc.line - lines_removed - 1, before + f"opp_par_loop_{loop.iterator_type.name}__{loop.kernel}" + after)

        elif loop.loop_type == OP.LoopType.MOVE_LOOP:
            before, after = buffer.get(loop.loc.line - lines_removed - 1).split("opp_particle_move", 1)
            after = re.sub(
                rf"{loop.kernel}\s*,\s*\"{loop.loop_name}\"\s*,\s*", "", after, count=1
            )  # TODO: This assumes that the kernel arg is on the same line as the call
            buffer.update(loop.loc.line - lines_removed - 1, before + f"opp_particle_move__{loop.kernel}" + after)

    # 3. Update headers
    index = buffer.search(r'\s*#include\s+"opp_templates\.h"')
    assert index is not None
    index += 2

    # buffer.insert(index, '#ifdef OPENACC\n#ifdef __cplusplus\nextern "C" {\n#endif\n#endif\n')
    loop_list = []
    move_loop = None
    for loop in program.loops:
        if loop.loop_type == OP.LoopType.PAR_LOOP:
            prototype = f'void opp_par_loop_{loop.iterator_type.name}__{loop.kernel}(opp_set,opp_iterate_type{",opp_arg" * len(loop.args)});'
        elif loop.loop_type == OP.LoopType.MOVE_LOOP:
            prototype = f'void opp_particle_move__{loop.kernel}(opp_set,opp_map,opp_map{",opp_arg" * len(loop.args)});'
            move_loop = loop

        if prototype not in loop_list:
            loop_list.append(prototype)
            buffer.insert(index, prototype)

    pattern = r"opp_init_direct_hop\([^)]*\)"
    dh_api_match = re.search(pattern, source)

    if dh_api_match:
        prototype = f'void opp_init_direct_hop_cg(double,int,const opp_dat,const opp::BoundingBox&,opp_map,opp_map{",opp_arg" * len(move_loop.args)});'
        buffer.insert(index, prototype)
        index += 1
        move_loop.dh_loop_required = True

    source = buffer.translate()

    # Substitude dh init dats and arguments
    if dh_api_match:
        source = re.sub(r'opp_init_direct_hop\(([^)]*)\);', \
                    r'opp_init_direct_hop_cg(\g<1>{});'.format(get_dh_init_args(move_loop)), source)
    
    # Substitude the opp_templates.h include for opp_lib.h
    source = re.sub(r'#include\s+"opp_templates\.h"', '#include "opp_lib.h"', source)

    source = re.sub(r'// USER WRITTEN CODE', '// AUTO GENERATED CODE', source)
    source = re.sub(r'#include "kernels.h"', '// #include "kenels.h" // codegen commented...', source)

    return source
