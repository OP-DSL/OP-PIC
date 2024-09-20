from __future__ import annotations

import dataclasses
from dataclasses import dataclass
from enum import Enum
from typing import TYPE_CHECKING, Callable, List, Optional, Union

from util import ABDC, findIdx

if TYPE_CHECKING:
    from store import Location


class AccessType(Enum):
    READ = 0
    WRITE = 1
    RW = 2

    INC = 3
    MIN = 4
    MAX = 5

    WORK = 6

    @staticmethod
    def values() -> List[int]:
        return [x.value for x in list(AccessType)]

    @staticmethod
    def names() -> List[str]:
        return [x.name for x in list(AccessType)]


class OpError(Exception):
    message: str
    loc: Optional[Location]

    def __init__(self, message: str, loc: Optional[Location] = None) -> None:
        self.message = message
        self.loc = loc

    def __str__(self) -> str:
        if self.loc is not None:
            return f"{self.loc}: OP error: {self.message}"
        else:
            return f"OP error: {self.message}"


class Type:
    formatter: Callable[["Type"], str]

    @classmethod
    def set_formatter(cls, formatter: Callable[["Type"], str]) -> None:
        cls.formatter = formatter

    def __str__(self) -> str:
        return self.__class__.formatter(self)

class IterateType(Enum):
    
    all = 1
    injected = 2

    @staticmethod
    def values() -> List[int]:
        return [x.value for x in list(IterateType)]

    @staticmethod
    def names() -> List[str]:
        return [x.name for x in list(IterateType)]
    
    def __str__(self) -> str:
        return self.name

class LoopType(Enum):
    
    PAR_LOOP = 1
    MOVE_LOOP = 2

    @staticmethod
    def values() -> List[int]:
        return [x.value for x in list(IterateType)]

    @staticmethod
    def names() -> List[str]:
        return [x.name for x in list(IterateType)]

class SetType(Enum):
    
    MESH = 1
    PARTICLE = 2

    @staticmethod
    def values() -> List[int]:
        return [x.value for x in list(SetType)]

    @staticmethod
    def names() -> List[str]:
        return [x.name for x in list(SetType)]
            
@dataclass(frozen=True)
class Int(Type):
    signed: bool
    size: int

    def __repr__(self) -> str:
        if self.signed and self.size == 32:
            return "int"
        elif self.size == 32:
            return "unsigned"
        else:
            return f"{'i' if self.signed else 'u'}{self.size}"

    def __str__(self) -> str:
        return "OPP_INT"

    # @classmethod
    # def int_formatter(cls, int_instance: "Int") -> str:
    #     return "OPP_INT"
        # if int_instance.signed and int_instance.size == 32:
        #     return "int"
        # elif int_instance.size == 32:
        #     return "unsigned"
        # else:
        #     return f"{'i' if int_instance.signed else 'u'}{int_instance.size}"
    
# Int.set_formatter(Int.int_formatter)

@dataclass(frozen=True)
class Float(Type):
    size: int

    def __repr__(self) -> str:
        if self.size == 32:
            return "float"
        elif self.size == 64:
            return "double"
        else:
            return "f{self.size}"

    def __str__(self) -> str:
        return "OPP_REAL"

    # @classmethod
    # def float_formatter(cls, float_instance: "Float") -> str:
    #     return "OPP_REAL"
        # if float_instance.size == 32:
        #     return "float"
        # elif float_instance.size == 64:
        #     return "double"
        # else:
        #     return "f{float_instance.size}"
    
# Float.set_formatter(Float.float_formatter)

@dataclass(frozen=True)
class Bool(Type):
    pass

    def __repr__(self) -> str:
        return "bool"

    def __str__(self) -> str:
        return "OPP_BOOL"

@dataclass(frozen=True)
class Custom(Type):
    name: str

    def __repr__(self) -> str:
        return self.name


@dataclass(frozen=True)
class Const:
    loc: Location
    ptr: str

    dim: int
    typ: Type

    def __str__(self) -> str:
        return f"Const(loc={self.loc}, ptr='{self.ptr}', dim={self.dim}, typ={self.typ})"

@dataclass
class OppSet:
    id: int
    loc: Location
    name: str
    set_type: SetType
    cell_set: Optional[OppSet]

    def __init__(self, id : int, loc: Location, name: str, set_type: SetType, cell_set: OppSet = None) -> None:
        self.id = id
        self.loc = loc
        self.name = name
        self.set_type = set_type
        self.cell_set = cell_set
    
    def __str__(self) -> str:
        return f"Set(id={self.id}, name='{self.name}', set_type={self.set_type}, cell_set={self.cell_set}, loc={self.loc})"

@dataclass
class Map:
    id: int

    ptr: str
    arg_id: int
    dim: Optional[int]
    from_set: Optional[OppSet]
    to_set: Optional[OppSet]
    loc: Optional[Location]

    def __init__(self, id: int, ptr: str, arg_id: int, dim: int, from_set: OppSet, to_set: OppSet, loc: Optional[Location]) -> None:
        self.id = id
        self.ptr = ptr
        self.arg_id = arg_id
        self.loc = loc
        self.dim = dim
        self.from_set = from_set
        self.to_set = to_set
    
    def __str__(self) -> str:
        return f"Map(id={self.id}, ptr='{self.ptr}', dim={self.dim}, from_set={self.from_set}, to_set={self.to_set}, loc={self.loc})"

@dataclass
class Dat:
    id: int

    ptr: str
    arg_id: int

    dim: Optional[int]
    typ: Type
    soa: bool

    set: Optional[OppSet]
    loc: Optional[Location]
    flag = True

    def __str__(self) -> str:
        return (
            f"Dat(id={self.id}, ptr='{self.ptr}', arg_id={self.arg_id}, dim={self.dim}, typ={self.typ}, soa={self.soa}, set={self.set}, loc={self.loc}, flag={self.flag}))"
        )


@dataclass
class Arg(ABDC):
    id: int
    loc: Location


@dataclass
class ArgDat(Arg):
    access_type: AccessType
    opt: bool

    dat_id: int

    map_id: Optional[int]
    map_idx: Optional[int]
    p2c_id: Optional[int]
    offset: bool
    flag = True

    def __str__(self) -> str:
        return (
            f"ArgDat(id={self.id}, loc={self.loc}, access_type={str(self.access_type) + ',':17} opt={self.opt}, "
            f"dat_id={self.dat_id}, map_id={self.map_id}, map_idx={self.map_idx}, p2c_id={self.p2c_id}, flag={self.flag}, offset={self.offset})"
        )


@dataclass
class ArgGbl(Arg):
    access_type: AccessType
    opt: bool

    ptr: str

    dim: Optional[int]
    typ: Type

    def __str__(self) -> str:
        return (
            f"ArgGbl(id={self.id}, loc={self.loc}, access_type={str(self.access_type) + ',':17} opt={self.opt}, "
            f"ptr={self.ptr}, dim={self.dim}, typ={self.typ})"
        )


@dataclass
class ArgIdx(Arg):
    map_id: Optional[int]
    map_idx: Optional[int]

    def __str__(self) -> str:
        return f"ArgIdx(id={self.id}, loc={self.loc}, map_id={self.map_id}, map_idx={self.map_idx})"


@dataclass
class ArgInfo(Arg):
    ptr: str

    dim: Optional[int]
    typ: Type

    ref: int

    def __str__(self) -> str:
        return f"ArgInfo(id={self.id}, loc={self.loc}, ptr={self.ptr}, dim={self.dim}, typ={self.typ}, ref={self.ref})"


class Loop:
    name: str

    loc: Location
    kernel: str

    iterator_set : OppSet
    iterator_type : IterateType
    loop_type : LoopType
    loop_name : str
    p2c_map : Optional[Dat]
    c2c_map : Optional[Map]

    args: List[Arg]
    args_expanded: List[Arg]

    dats: List[Dat]
    maps: List[Map]

    consts: Set[Const]

    fallback: bool
    dh_loop_required = False
    
    def __init__(self, name: str, loc: Location, kernel: str, iterator_set : OppSet, iterator_type : IterateType, loop_name: str, loop_type = LoopType.PAR_LOOP, p2c_map=None, c2c_map = None) -> None:
        self.name = name

        self.loc = loc
        self.kernel = kernel
        self.iterator_set = iterator_set
        self.iterator_type = iterator_type
        self.loop_type = loop_type
        self.loop_name = loop_name
        self.p2c_map = p2c_map
        self.c2c_map = c2c_map

        self.dats = []
        self.maps = []

        self.args = []
        self.args_expanded = []

        self.consts = set()

        self.fallback = False

    def addArgDat(
        self,
        loc: Location,
        dat_ptr: str,
        map_ptr: Optional[str],
        map_idx: Optional[int],
        p2c_ptr: Optional[str],
        access_type: AccessType,
        offset: bool,
        opt: bool = False,
    ) -> None:
        arg_id = len(self.args)

        dat_dim = 22
        dat_typ = 22
        dat_soa = 22

        dat_id = findIdx(self.dats, lambda d: d.ptr == dat_ptr)
        if dat_id is None:
            dat_id = len(self.dats)

            dat = Dat(dat_id, dat_ptr, arg_id, dat_dim, dat_typ, dat_soa, None, None)
            self.dats.append(dat)
        elif self.dats[dat_id].dim is None and dat_dim is not None:
            self.dats[dat_id] = dataclasses.replace(self.dats[dat_id], dim=dat_dim)

        # p2c_id = p2c_ptr
        # if p2c_ptr is not None:
        #     p2c_id = findIdx(self.dats, lambda d: d.ptr == p2c_ptr)
        #     if p2c_id is None:
        #         p2c_id = len(self.dats)

        #         p2c_dat = Dat(p2c_id, p2c_ptr, arg_id, dat_dim, dat_typ, dat_soa, None, None)
        #         self.dats.append(p2c_dat)
        #     elif self.dats[p2c_id].dim is None and dat_dim is not None:
        #         self.dats[p2c_id] = dataclasses.replace(self.dats[p2c_id], dim=dat_dim)

        map_id = None
        if map_ptr is not None:
            map_id = findIdx(self.maps, lambda m: m.ptr == map_ptr)

            if map_id is None:
                map_id = len(self.maps)

                map_ = Map(map_id, map_ptr, arg_id, None, None, None, None)
                self.maps.append(map_)

        arg = ArgDat(arg_id, loc, access_type, opt, dat_id, map_id, map_idx, p2c_ptr, offset)
        self.args.append(arg)

        if map_ptr is None or map_idx is None or map_idx >= 0:
            self.args_expanded.append(arg)
            return

        for real_map_idx in range(-map_idx):
            arg_expanded = dataclasses.replace(arg, map_idx=real_map_idx)
            self.args_expanded.append(arg_expanded)

    def addArgGbl(
        self, loc: Location, ptr: str, dim: Optional[int], typ: Type, access_type: AccessType, opt: bool
    ) -> None:
        arg_id = len(self.args)
        arg = ArgGbl(arg_id, loc, access_type, opt, ptr, dim, typ)

        self.args.append(arg)
        self.args_expanded.append(arg)

    def addArgIdx(self, loc: Location, map_ptr: Optional[str], map_idx: Optional[int]) -> None:
        arg_id = len(self.args)

        map_id = None
        if map_ptr is not None:
            map_id = findIdx(self.maps, lambda m: m.ptr == map_ptr)

            if map_id is None:
                map_id = len(self.maps)

                map_ = Map(map_id, map_ptr, arg_id, None, None, None, None)
                self.maps.append(map_)

        arg = ArgIdx(arg_id, loc, map_id, map_idx)
        self.args.append(arg)

    def addArgInfo(self, loc: Location, ptr: str, dim: Optional[int], typ: Type, ref: int) -> None:
        arg_id = len(self.args)
        arg = ArgInfo(arg_id, loc, ptr, dim, typ, ref)

        self.args.append(arg)
        self.args_expanded.append(arg)

    def optIdx(self, arg: Arg) -> Optional[int]:
        idx = 0
        for arg2 in self.args:
            if arg2 == arg:
                break

            if getattr(arg2, "opt", None) is not None:
                idx += 1

        return idx

    def addConst(self, const: Const) -> None:
        self.consts.add(const)

    def arg(self, x: Union[Dat, int]) -> Optional[Arg]:
        if isinstance(x, Dat):
            return self.args[x.arg_id]

        if isinstance(x, int) and x < len(self.args):
            return self.args[x]

        return None

    def dat(self, x: Union[ArgDat, int]) -> Optional[Dat]:
        if isinstance(x, Arg):
            return self.dats[x.dat_id]

        if isinstance(x, int) and x < len(self.dats):
            return self.dats[x]

        return None

    def map(self, x: Union[ArgDat, int]) -> Optional[Map]:
        if isinstance(x, ArgDat) and x.map_id is not None:
            return self.maps[x.map_id]

        if isinstance(x, int) and x < len(self.maps):
            return self.maps[x]

        return None

    def __str__(self) -> str:
        args = "\n    ".join([str(a) for a in self.args])
        argsEx = "\n    ".join([str(a) for a in self.args_expanded])
        dat_str = "\n    ".join([str(d) for d in self.dats])
        map_str = "\n    ".join([str(m) for m in self.maps])
        consts_str = "\n    ".join([str(m) for m in self.consts])

        if len(self.dats) > 0:
            dat_str = f"\n    {dat_str}\n"

        if len(self.maps) > 0:
            map_str = f"\n    {map_str}\n"

        if len(self.consts) > 0:
            consts_str = f"\n    {consts_str}\n"

        return (
            f"Loop at {self.loc}:\n    Name: {self.name}\n    Kernel function: {self.kernel}\n    Iterating Set: {self.iterator_set}\n    p2c_map: {self.p2c_map}\n    c2c_map: {self.c2c_map} \n\n    args:\n    {args}\n    argsEx:\n    {argsEx}\n"
            + dat_str
            + map_str
            + consts_str
        )
