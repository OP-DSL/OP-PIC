import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

from clang.cindex import Cursor, CursorKind, TranslationUnit, TypeKind
from clang.cindex import conf as clang_internal  # type: ignore

import opp as OPP
from store import Function, Location, ParseError, Program, Type
from util import safeFind

def parseSet(node: Cursor, name: str, program: Program, set_type : OPP.SetType) -> None:

    children = list(node.get_children())
    if len(children) < 3:
        raise ParseError(f'unable to parse Set, argument mismatch')
        
    loc = parseLocation(node) # TODO : check whether the location is correct! not required for now
    cell_set = None
    id = len(program.sets)

    if (set_type == OPP.SetType.PARTICLE):
        cell_set = parseIdentifier(children[2])

    set = OPP.OppSet(id, loc, name, set_type, cell_set)
    program.sets.append(set)

    # print(f'parseSet {name} -> {set}')

def parseMap(node: Cursor, name: str, program: Program) -> None:

    children = list(node.get_children())
    if len(children) < 6:
        raise ParseError(f'unable to parse Map, argument mismatch')
    
    loc = parseLocation(node) # TODO : check whether the location is correct! not required for now
    from_set = parseIdentifier(children[1])
    to_set = parseIdentifier(children[2])
    dim = parseIntExpression(children[3])
    id = len(program.maps)

    map = OPP.Map(id, name, None, dim, from_set, to_set, loc)
    program.maps.append(map)

    # print(f'parseMap {name} -> {map}')

def parseDat(node: Cursor, name: str, program: Program) -> None:

    children = list(node.get_children())
    if len(children) < 6:
        raise ParseError(f'unable to parse Dat, argument mismatch')
    
    loc = parseLocation(node) # TODO : check whether the location is correct! not required for now
    set = parseIdentifier(children[1])
    dim = parseIntExpression(children[2])
    typ = parseIdentifier(children[3])

    if (str(typ) == "DT_INT"):
        typ = "OPP_INT" 
    elif (str(typ) == "DT_REAL"):
        typ = "OPP_REAL" 
    elif (str(typ) == "DT_BOOL"):
        typ = "OPP_BOOL" 
    else:
        raise ParseError(f'unable to parse Dat, unknown data type {typ}')
    
    id = len(program.dats)

    dat = OPP.Dat(id, name, None, dim, typ, True, set, loc)
    program.dats.append(dat)

    # print(f'parseDat {name} -> {dat}')

def parseMeta(node: Cursor, program: Program) -> None:
    if node.kind == CursorKind.TYPE_REF:
        parseTypeRef(node, program)

    if node.kind == CursorKind.FUNCTION_DECL:
        parseFunction(node, program)

    if node.kind == CursorKind.VAR_DECL:
        
        args = list(node.get_children())
        if len(args) < 2 or args[1].kind != CursorKind.CALL_EXPR:
            return

        children = list(args[1].get_children())
        if len(children) < 3:
            return
        
        child = list(children[0].get_children())
        if len(child) < 1:
            return

        if child[0].spelling == "opp_decl_set" or child[0].spelling == "opp_decl_set_hdf5":
            parseSet(args[1], node.spelling, program, OPP.SetType.MESH)
        elif child[0].spelling == "opp_decl_particle_set" or child[0].spelling == "opp_decl_particle_set_hdf5":
            parseSet(args[1], node.spelling, program, OPP.SetType.PARTICLE)
        elif child[0].spelling == "opp_decl_map" or child[0].spelling == "opp_decl_map_hdf5":
            parseMap(args[1], node.spelling, program)
        elif child[0].spelling == "opp_decl_dat" or child[0].spelling == "opp_decl_dat_hdf5":
            parseDat(args[1], node.spelling, program)

    for child in node.get_children():
        parseMeta(child, program)


def parseTypeRef(ref_node: Cursor, program: Program) -> None:
    node = ref_node.get_definition()

    if node is None or Path(str(node.location.file)) != program.path:
        return

    # print(f'ZAMO parseTypeRef {node.spelling}')
    
    matching_entities = program.findEntities(node.spelling)
    for entity in matching_entities:
        if entity.ast == node:
            return

    typ = Type(node.spelling, node, program)

    for n in node.walk_preorder():
        if n.kind != CursorKind.CALL_EXPR and n.kind != CursorKind.TYPE_REF:
            continue

        n = n.get_definition()

        if n is None or Path(str(n.location.file)) != program.path:
            continue

        typ.depends.add(n.spelling)

    program.entities.append(typ)


def parseFunction(ref_node: Cursor, program: Program) -> None:
    node = ref_node.get_definition()

    if node is None or Path(str(node.location.file)) != program.path:
        return

    matching_entities = program.findEntities(node.spelling)
    for entity in matching_entities:
        if entity.ast == node:
            return

    function = Function(node.spelling, node, program)

    for n in node.get_children():
        if n.kind != CursorKind.PARM_DECL:
            continue

        # param_type = n.type.get_canonical()
        # typ, _ = parseType(param_type.spelling, parseLocation(n), True)

        function.parameters.append(n.spelling)

    for n in node.walk_preorder():
        if n.kind != CursorKind.CALL_EXPR and n.kind != CursorKind.TYPE_REF:
            continue

        n = n.get_definition()

        if n is None or Path(str(n.location.file)) != program.path:
            continue

        function.depends.add(n.spelling)

    program.entities.append(function)


def parseLoops(translation_unit: TranslationUnit, program: Program) -> None:
    macros: Dict[Location, str] = {}
    nodes: List[Cursor] = []

    for node in translation_unit.cursor.get_children():
        if node.kind == CursorKind.MACRO_DEFINITION:
            continue
        elif node.location.file.name != translation_unit.spelling:  # type: ignore
            continue
        elif node.kind == CursorKind.MACRO_INSTANTIATION:
            macros[parseLocation(node)] = node.spelling
            continue
        
        # print(f"ZAM parseLoops appending nodes {node.spelling}")
        nodes.append(node)

    for node in nodes:
        for child in node.walk_preorder():
            if child.kind.is_unexposed():
            # if child.kind == CursorKind.CALL_EXPR:
                # print(f"ZAM parseLoops - calling parseCall node:{node.spelling} child:{child.spelling}")
                parseCall(child, macros, program)

def parseUnexposedFunction(node: Cursor) -> Union[Tuple[str, List[Cursor]], None]:
    args = []
    for child in node.get_children():
        args.append(child)

    if len(args) == 0:
        return None

    first_child = args.pop(0)

    if (first_child.kind == CursorKind.MEMBER_REF_EXPR and len(list(first_child.get_children())) >= 2):     
        name_token = list(first_child.get_children())[1]
        name = name_token.spelling
    elif first_child.kind == CursorKind.DECL_REF_EXPR and len(list(first_child.get_tokens())) > 0:
        name = list(first_child.get_tokens())[0].spelling
    else:
        return None

    return (name, args)

def parseCall(node: Cursor, macros: Dict[Location, str], program: Program) -> None:

    if parseUnexposedFunction(node) == None:
        return
    else:
        (name, args) = parseUnexposedFunction(node)

    loc = parseLocation(node)

    if name == "opp_decl_const":
        program.consts.append(parseConst(args, loc))
    elif name == "opp_par_loop":
        program.loops.append(parseLoop(program, args, loc, macros))
    elif name == "opp_particle_move":
        program.loops.append(parseMoveLoop(program, args, loc, macros))

def parseConst(args: List[Cursor], loc: Location) -> OPP.Const:
    if len(args) != 3:
        raise ParseError("incorrect number of args passed to op_decl_const", loc)

    # TODO dim may not be literal
    dim = parseIntExpression(args[0])
    name = args[2].spelling[1:-1]
    
    typ_str = str()
    if "double" in args[1].type.spelling:
        typ_str = "double"
    elif "int" in args[1].type.spelling:
        typ_str = "int"
    elif "bool" in args[1].type.spelling:
        typ_str = "bool"
    else:
        raise ParseError("incorrect type passed to opp_decl_const", loc)
        
    typ, _ = parseType(typ_str, loc)

    # print(f"parseConst {dim} {name} {typ}")

    return OPP.Const(loc, name, dim, typ)

def parseLoop(program: Program, args: List[Cursor], loc: Location, macros: Dict[Location, str]) -> OPP.Loop:
    if len(args) < 4:
        raise ParseError("incorrect number of args passed to opp_par_loop")

    kernel = parseIdentifier(args[0])
    loop_name = args[1].spelling[1:-1]
    iter_set = parseIdentifier(args[2])
    name = f"{program.path.stem}_{len(program.loops) + 1}_{kernel}"
    iterate_type = parseIterateType(args[3])

    loop = OPP.Loop(name, loc, kernel, iter_set, iterate_type, loop_name)

    for node in args[4:]:
        name = node.spelling

        arg_loc = parseLocation(node)
        arg_args = list(node.get_arguments())

        # print(f'ZAM node.spelling {node.spelling} node.kind {node.kind} arg_loc {arg_loc}')
        # for arg in arg_args:
        #     print(f'ZAM opp_arg_dat args -> name {arg.spelling} kind {arg.kind}')

        if name == "op_arg_idx":
            parseArgIdx(loop, arg_args, arg_loc, macros)
            
        elif name == "opp_arg_dat":
            parseArgDat(loop, False, arg_args, arg_loc, macros)

        elif name == "opp_arg_gbl":
            parseArgGbl(loop, False, arg_args, arg_loc, macros)

        else:
            raise ParseError(f"invalid loop argument {name}", parseLocation(node))

    return loop

def parseMoveLoop(program: Program, args: List[Cursor], loc: Location, macros: Dict[Location, str]) -> OPP.Loop:
    if len(args) < 5:
        raise ParseError("incorrect number of args passed to opp_particle_move")

    kernel = parseIdentifier(args[0])
    loop_name = args[1].spelling[1:-1]
    iter_set = parseIdentifier(args[2])
    c2c_map = parseIdentifier(args[3])
    p2c_map = parseIdentifier(args[4])
    name = f"{program.path.stem}_{len(program.loops) + 1}_{kernel}"
    
    loop = OPP.Loop(name, loc, kernel, iter_set, OPP.IterateType.all, loop_name, OPP.LoopType.MOVE_LOOP, p2c_map, c2c_map)

    for node in args[5:]:
        # node = descend(descend(node))
        name = node.spelling

        arg_loc = parseLocation(node)
        arg_args = list(node.get_arguments())

        # print(f'ZAM parseMoveLoop node.spelling {node.spelling} node.kind {node.kind} arg_loc {arg_loc}')
        # for arg in arg_args:
        #     print(f'ZAM parseMoveLoop args -> name {arg.spelling} kind {arg.kind}')

        if name == "op_arg_idx":
            parseArgIdx(loop, arg_args, arg_loc, macros)
            
        elif name == "opp_arg_dat":
            parseArgDat(loop, False, arg_args, arg_loc, macros)

        elif name == "opp_arg_gbl":
            parseArgGbl(loop, False, arg_args, arg_loc, macros)

        else:
            raise ParseError(f"invalid loop argument {name}", parseLocation(node))

    return loop

def parseArgDat(loop: OPP.Loop, opt: bool, args: List[Cursor], loc: Location, macros: Dict[Location, str]) -> None:
    
    if len(args) < 2 or len(args) > 6:
        raise ParseError("incorrect number of args passed to opp_arg_dat", loc)

    dat_ptr = parseIdentifier(args[0])
    access_type = parseAccessType(args[len(args)-2], loc, macros)
    
    offset = True
    if args[len(args)-1].kind == CursorKind.CXX_BOOL_LITERAL_EXPR:
        offset = parseBoolExpression(args[len(args)-1])

    map_idx = -1
    map_ptr = None
    p2c_ptr = None

    if len(args) == 4: # map from particle to cell index     
        p2c_ptr = parseIdentifier(args[1])
    elif len(args) == 5: # indirect mapping     
        map_idx = parseIntExpression(args[1])
        map_ptr = parseIdentifier(args[2])
    elif len(args) == 6: # double indirect mappings
        map_idx = parseIntExpression(args[1])
        map_ptr = parseIdentifier(args[2])
        p2c_ptr = parseIdentifier(args[3])

    # print(f'ZAM parseArgDat Len {len(args)} | {loc} | {dat_ptr} {map_idx} {map_ptr} {p2c_ptr} {access_type} | {offset} | {args[len(args)-1].kind}')

    loop.addArgDat(loc, dat_ptr, map_ptr, map_idx, p2c_ptr, access_type, offset)


def parseArgGbl(loop: OPP.Loop, opt: bool, args: List[Cursor], loc: Location, macros: Dict[Location, str]) -> None:
    if len(args) != 4:
        raise ParseError("incorrect number of args passed to opp_arg_gbl", loc)

    ptr = parseIdentifier(args[0])
    dim = parseIntExpression(args[1])
    typ, _ = parseType(parseStringLit(args[2]), loc)

    access_type = parseAccessType(args[3], loc, macros)

    # print(f'ZAM parseArgGbl Len {len(args)} | {loc} | {ptr} {dim} {access_type}')

    loop.addArgGbl(loc, ptr, dim, typ, access_type, opt)


def parseArgIdx(loop: OPP.Loop,  args: List[Cursor], loc: Location, macros: Dict[Location, str]) -> None:
    if args is None or len(args) != 2:
        raise ParseError("incorrect number of arguments for op_arg_idx", loc)

    map_idx = parseIntExpression(args[0])
    map_ptr = None if macros.get(parseLocation(args[1])) == "OP_ID" else parseIdentifier(args[1])

    loop.addArgIdx(loc, map_ptr, map_idx)


def parseIdentifier(node: Cursor, raw: bool = True) -> str:
    if raw:
        return "".join([t.spelling for t in node.get_tokens()])

    while node.kind == CursorKind.CSTYLE_CAST_EXPR:
        node = list(node.get_children())[1]

    if node.kind == CursorKind.UNEXPOSED_EXPR:
        node = descend(node)

    if node.kind == CursorKind.UNARY_OPERATOR and next(node.get_tokens()).spelling in ("&", "*",):
        node = descend(node)

    if node.kind == CursorKind.GNU_NULL_EXPR:
        raise ParseError("expected identifier, found NULL", parseLocation(node))

    if node.kind != CursorKind.DECL_REF_EXPR:
        raise ParseError("expected identifier", parseLocation(node))

    return node.spelling


def parseIntExpression(node: Cursor) -> int:
    if not (node.type.kind != TypeKind.INT or node.type.kind != TypeKind.ENUM):
        raise ParseError("expected int expression", parseLocation(node))

    eval_result = clang_internal.lib.clang_Cursor_Evaluate(node)
    val = clang_internal.lib.clang_EvalResult_getAsInt(eval_result)
    clang_internal.lib.clang_EvalResult_dispose(eval_result)

    return val

def parseStringLit(node: Cursor) -> str:
    if node.kind != CursorKind.UNEXPOSED_EXPR:
        raise ParseError("expected string literal")

    node = descend(node)
    if node.kind != CursorKind.STRING_LITERAL:
        raise ParseError("expected string literal")

    return node.spelling[1:-1]

def parseBoolExpression(node: Cursor) -> bool:
    if node.kind != CursorKind.CXX_BOOL_LITERAL_EXPR:
        raise ParseError("expected boolean expression", parseLocation(node))

    eval_result = clang_internal.lib.clang_Cursor_Evaluate(node)
    val = clang_internal.lib.clang_EvalResult_getAsInt(eval_result) != 0
    clang_internal.lib.clang_EvalResult_dispose(eval_result)

    return val

def parseAccessType(node: Cursor, loc: Location, macros: Dict[Location, str]) -> OPP.AccessType:
    access_type_raw = parseIntExpression(node)

    if access_type_raw not in OPP.AccessType.values():
        raise ParseError(
            f"invalid access type {access_type_raw}, expected one of {', '.join(OPP.AccessType.names())}", loc)

    return OPP.AccessType(access_type_raw)

def parseIterateType(node: Cursor) -> OPP.IterateType:
    iteration_type_raw = parseIntExpression(node)

    if iteration_type_raw not in OPP.IterateType.values():
        raise ParseError(
            f"invalid iterator type {iteration_type_raw}, expected one of {', '.join(OPP.IterateType.names())}")

    return OPP.IterateType(iteration_type_raw)

def parseType(typ: str, loc: Location, include_custom=False) -> Tuple[OPP.Type, bool]:
    typ_clean = typ.strip()
    typ_clean = re.sub(r"\s*const\s*", "", typ_clean)

    soa = False
    if re.search(r":soa", typ_clean):
        soa = True

    typ_clean = re.sub(r"\s*:soa\s*", "", typ_clean)

    typ_map = {
        "int": OPP.Int(True, 32),
        "uint": OPP.Int(False, 32),
        "ll": OPP.Int(True, 64),
        "ull": OPP.Int(False, 64),
        "float": OPP.Float(32),
        "double": OPP.Float(64),
        "bool": OPP.Bool(),
    }

    if typ_clean in typ_map:
        return typ_map[typ_clean], soa

    if include_custom:
        return OPP.Custom(typ_clean), soa

    raise ParseError(f'unable to parse type "{typ}"', loc)


def parseLocation(node: Cursor) -> Location:
    return Location(node.location.file.name, node.location.line, node.location.column)  # type: ignore


def descend(node: Cursor) -> Cursor:
    return next(node.get_children())


def descend_opt(node: Optional[Cursor]) -> Optional[Cursor]:
    if node is None:
        return None

    return next(node.get_children(), None)