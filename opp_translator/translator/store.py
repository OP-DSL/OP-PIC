from __future__ import annotations

import os
import copy
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, List, Optional, Set, Tuple
from util import flatten, uniqueBy
from clang.cindex import Cursor, CursorKind, SourceRange

# import ops
import op as OP
from op import OpError as OpsError
from util import ABDC, findIdx

if TYPE_CHECKING:
    from language import Lang


@dataclass(frozen=True)
class Location:
    file: str
    line: int
    column: int

    def __str__(self) -> str:
        return f"{os.path.abspath(self.file)}/{self.line}:{self.column}"


@dataclass
class ParseError(Exception):
    message: str
    loc: Optional[Location] = None

    def __str__(self) -> str:
        if self.loc:
            return f"Parse error at {self.loc}: {self.message}"
        else:
            return f"Parse error: {self.message}"

@dataclass
class Entity:
    name: str
    ast: Any

    program: Program
    scope: List[str] = field(default_factory=list)
    depends: Set[str] = field(default_factory=set)

    #Deep-copy everything but the program reference
    def __deepcopy__(self, memo) -> Optional["Entity"]:
        cls = self.__class__
        result = cls.__new__(cls)

        setattr(result, "program", self.program)

        memo[id(self)] = result

        for k, v in self.__dict__.items():
            if hasattr(result, k):
                continue

            setattr(result, k, copy.deepcopy(v, memo))

        return result

@dataclass
class Type(Entity):
    def __str__(self) -> str:
        return f"Type(name='{self.name}', scope={self.scope}, depends={self.depends})"

@dataclass
class Function(Entity):
    parameters: List(str) = field(default_factory=list)
    returns: Optional[OP.Type] = None

    def __str__(self) -> str:
        return f"Function(name='{self.name}', scope={self.scope}, depends={self.depends})"

@dataclass
class Program:
    path: Path

    ast: Any
    source: str

    consts: List[OP.Const] = field(default_factory=list)
    loops: List[OP.Loop] = field(default_factory=list)

    entities: List[Entity] = field(default_factory=list)

    sets: List[OP.Set] = field(default_factory=list)
    maps: List[OP.Map] = field(default_factory=list)
    dats: List[OP.Dat] = field(default_factory=list)

    ndim: Optional[int] = None
    soa_val: Optional[bool] = False
    init_flag: Optional[bool] = False

    def findEntities(self, name: str, scope: List[str] = []) -> List[Entity]:
        def in_scope(entity: Entity):
            return len(entity.scope) <= len(scope) and all(map(lambda s1, s2: s1 == s2, zip(entity.scope, scope)))

        candidates = list(filter(lambda e: e.name == name and in_scope(e), self.entities))

        if len(candidates) == 0:
            return []

        candidates.sort(key=lambda e: len(e.scope), reverse=True)
        min_scope = len(candidates[0].scope)

        #returning canditages with min scope    
        return list(filter(lambda e: len(e.scope) == min_scope, candidates))

    def __str__(self) -> str:
        outString = "\nprogram path=" + str(self.path)  + ",\n"
        outString += "ast=" + str(self.ast) + ",\n"
        outString += "ndim=" + str(self.ndim) + ",\n"

        outString += "\n---------------------\n"
        outString += "       consts        \n"
        outString += "---------------------\n"
        for const in self.consts:
            outString += str(const) + "\n"

        outString += "\n---------------------\n"    
        outString += "        loops        \n"
        outString += "---------------------\n"
        for loop in self.loops:
            outString += str(loop) + "\n"

        outString += "\n---------------------\n"    
        outString += "       Entities      \n"
        outString += "---------------------\n"
        for entity in self.entities:
            outString += str(entity) + "\n"  
        return outString

    def enrichLoopData(self) -> None:
        # for loop in self.loops:
        #     print(loop)
        # for loop in self.sets:
        #     print(loop) 
        # for loop in self.dats:
        #     print(loop)          
        # for loop in self.maps:
        #     print(loop) 

        # updating the loop with the consts accessed within the kernel
        for loop in self.loops:
            entities = self.findEntities(loop.kernel, [])
            # TODO : Use depencies as well
            for entity in entities:
                for node in entity.ast.walk_preorder():
                    if node.kind != CursorKind.DECL_REF_EXPR:
                        continue
                    const_id = findIdx(self.consts, lambda d: d.ptr == node.spelling)
                    if const_id is not None:
                        loop.addConst(self.consts[const_id])
        for loop in self.loops:
            for arg in loop.args:
                if isinstance(arg, OP.ArgDat):
                    stored_dat_id = findIdx(self.dats, lambda m: m.ptr == loop.dats[arg.dat_id].ptr) 
                    if stored_dat_id is None:
                        stored_map_id = findIdx(self.maps, lambda m: m.ptr == loop.dats[arg.dat_id].ptr) 
                        if stored_map_id is not None:
                            arg.flag = False
                        else:
                            raise OpsError(f"enrichLoopData : {arg} is not a dat nor map")
        for loop in self.loops:
            for map in loop.maps:
                map_stored = self.maps[findIdx(self.maps, lambda m: m.ptr == map.ptr)]         
                map.dim = map_stored.dim
                map.from_set = map_stored.from_set
                map.to_set = map_stored.to_set
                map.loc = map_stored.loc

            for obj in loop.dats: 
                dat_idx = findIdx(self.dats, lambda d: d.ptr == obj.ptr)
                if dat_idx is not None:
                    dat_stored = self.dats[dat_idx]           
                    obj.dim = dat_stored.dim
                    obj.typ = dat_stored.typ
                    obj.soa = dat_stored.soa
                    obj.set = dat_stored.set
                    obj.loc = dat_stored.loc
                else:
                    map_idx = findIdx(self.maps, lambda d: d.ptr == obj.ptr)
                    if map_idx is not None:
                        map_stored = self.maps[map_idx]           
                        obj.dim = map_stored.dim
                        obj.typ = "OPP_INT"
                        obj.soa = True
                        obj.set = map_stored.from_set
                        obj.loc = map_stored.loc    
                        obj.flag = False                 
                    else:
                        raise OpsError(f"enrichLoopData : {obj} is not a dat nor map")
            loop.iterator_set = self.sets[findIdx(self.sets, lambda m: m.name == loop.iterator_set)]
            if loop.p2c_map not in [None, "nullptr", "NULL"] :
                loop.p2c_map = self.maps[findIdx(self.maps, lambda m: m.ptr == loop.p2c_map)]
            elif loop.loop_type == OP.LoopType.MOVE_LOOP:
                raise OpsError(f"enrichLoopData : {loop.name} require a valid p2c map")
            if loop.c2c_map not in [None, "nullptr", "NULL"] :
                loop.c2c_map = self.maps[findIdx(self.maps, lambda m: m.ptr == loop.c2c_map)]
            elif loop.loop_type == OP.LoopType.MOVE_LOOP:
                raise OpsError(f"enrichLoopData : {loop.name} require a valid c2c map")

@dataclass
class Application:
    programs: List[Program] = field(default_factory=list)
    global_dim: Optional[int] = None
    external_consts: Set[str] = field(default_factory=set)
    
    def __str__(self) -> str:
        if len(self.programs) > 0:
            programs_str = "\n".join([str(p) for p in self.programs])
        else:
            programs_str = "No programs"

        return programs_str

    def findEntities(self, name: str, program: Program = None, scope: List[str] = []) -> List[Entity]:
        candidates = []

        if program is not None:
            candidates = program.findEntities(name, scope)

        if len(candidates) > 0:
             return candidates

        for program2 in self.programs:
            if program2 == program:
                continue

            candidates = program2.findEntities(name)
            if len(candidates) > 0:
                break

        return candidates

    # All const ptrs including external_consts
    def constPtrs(self) -> Set[str]:
        ptrs = set(self.external_consts)

        for program in self.programs:
            ptrs.update(const.ptr for const in program.consts)

        return ptrs
    
    def consts(self) -> List[OP.Const]:
        consts = flatten(program.consts for program in self.programs)
        return uniqueBy(consts, lambda c: c.ptr)

    def loops(self) -> List[Tuple[OP.Loop, Program]]:
        return flatten(map(lambda l: (l, p), p.loops) for p in self.programs)

    def uniqueLoops(self) -> List[OP.Loop]:
        return uniqueBy(self.loops(), lambda m: m[0].kernel)
        for p in self.programs:
            id = findId

    def validate(self, lang: Lang) -> None:
        self.validateConst(lang)
        self.validateLoops(lang)

        for program in self.programs:
            if self.global_dim == None:
                self.global_dim = program.ndim
            elif program.ndim == None:
                program.ndim = self.global_dim
            elif self.global_dim != program.ndim:
                raise OpsError(f"ndim mismatch with global dim={self.global_dim} and program dim={program.ndim} of program={program.path}")


    def validateConst(self, lang: Lang) -> None:
        seen_const_ptrs: Set[str] = set()

        for program in self.programs:
            for const in program.consts:
                if const.ptr in seen_const_ptrs:
                    raise OpsError(f"Duplicate const declaration: {const.ptr}", const.loc)

                seen_const_ptrs.add(const.ptr)

                # if const.dim.isdigit() and int(const.dim) < 1:
                #     raise OpsError(f"Invalid const dimension: {const.dim} of const: {const.ptr}", const.loc)


    def validateLoops(self, lang: Lang) -> None:
        for loop, Program in self.loops():
            num_opts = len([arg for arg in loop.args if getattr(arg, "opt", False)])
            #if num_opts > 128:
            #    raise OpsError(f"number of optional arguments exceeds 128: {num_opts}", loop.loc)
            for arg in loop.args:
                if isinstance(arg, OP.ArgDat):
                    self.validateArgDat(arg, loop, lang)

                if isinstance(arg, OP.ArgGbl):
                    self.validateArgGbl(arg, loop, lang)

            #self.validateKernel(loop, program, lang) TODO


    def validateArgDat(self, arg: OP.ArgDat, loop: OP.Loop, lang: Lang) -> None:
        valid_access_types = [
            OP.AccessType.READ, 
            OP.AccessType.WRITE, 
            OP.AccessType.RW, 
            OP.AccessType.INC
            ]

        if arg.access_type not in valid_access_types:
            raise OpsError(f"Invalid access type for dat argument: {arg.access_type}, arg: {arg}", arg.loc)


    def validateArgGbl(self, arg: OP.ArgGbl, loop: OP.Loop, lang: Lang) -> None:
        valid_access_types = [
            OP.AccessType.READ, 
            OP.AccessType.WRITE, 
            OP.AccessType.RW, 
            OP.AccessType.INC,
            OP.AccessType.MIN,
            OP.AccessType.MAX
            ]

        if arg.access_type not in valid_access_types:
            raise OpsError(f"Invalid access type for gbl argument: {arg.access_type}", arg.loc)

        if arg.access_type in [OP.AccessType.INC, OP.AccessType.MIN, OP.AccessType.MAX] and arg.typ not in \
            [OP.Float(64), OP.Float(32), OP.Float(16), OP.Int(True, 32), OP.Int(False, 32), OP.Bool]:
            raise OpsError(f"Invalid access type for reduced gbl argument: {arg.access_type}", arg.loc)

        if str(arg.dim).isdigit() and int(str(arg.dim)) < 1:
            raise OpsError(f"Invalid gbl argument dimension: {arg.dim}", arg.loc)
 





