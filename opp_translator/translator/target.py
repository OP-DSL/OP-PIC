from __future__ import annotations

from typing import Any, Dict
from util import Findable

class Target(Findable):
    name: str
    kernel_translation: bool
    config: Dict[str, Any]

    def defaultConfig(self) -> Dict[str, Any]:
        return {}
    
    def __str__(self) -> str:
        return self.name

    def __eq__(self, other) -> bool:
        return self.name == other.name if type(other) is type(self) else False

    def __hash__(self) -> int:
        return hash(self.name)

    def matches(self, key: str) -> bool:
        return self.name == key.lower()

class Seq(Target):
    name = "seq"
    suffix = "seq"
    kernel_translation = False
    config = {"grouped": False, "device": 1, "omp": False}

    def defaultConfig(self) -> Dict[str, Any]:
        return {"grouped": False, "device": 1, "omp": False}

class Mpi(Target):
    name = "mpi"
    suffix = "mpi"
    kernel_translation = False
    config = {"grouped": False, "device": 1, "omp": False}

    def defaultConfig(self) -> Dict[str, Any]:
        return {"grouped": False, "device": 1, "omp": False}

class Omp(Target):
    name = "omp"
    suffix = "omp"
    kernel_translation = False
    config = {"grouped": False, "device": 1, "omp": True}

    def defaultConfig(self) -> Dict[str, Any]:
        return {"grouped": False, "device": 1, "omp": True}

class Cuda(Target):
    name = "cuda"
    suffix = "cuda"
    kernel_translation = True
    config = {"grouped": True, "device": 2, "atomics": True, "seg_red": False, "color2": False, "gbl_inc_atomic": False, "omp": False}
    
    def defaultConfig(self) -> Dict[str, Any]:
        return {"grouped": True, "device": 2, "atomics": True, "seg_red": False, "color2": False, "gbl_inc_atomic": False, "omp": False}

class Hip(Target):
    name = "hip"
    suffix = "hip"
    kernel_translation = True
    config = {"grouped": True, "device": 2, "atomics": True, "seg_red": True, "color2": False, "gbl_inc_atomic": False, "omp": False}
    
    def defaultConfig(self) -> Dict[str, Any]:
        return {"grouped": True, "device": 2, "atomics": True, "seg_red": True, "color2": False, "gbl_inc_atomic": False, "omp": False}

class Sycl(Target):
    name = "sycl"
    suffix = "sycl"
    kernel_translation = True
    config = {"grouped": True, "device": 2, "atomics": True, "color2": False, "gbl_inc_atomic": False, "omp": False}
    
    def defaultConfig(self) -> Dict[str, Any]:
        return {"grouped": True, "device": 2, "atomics": True, "color2": False, "gbl_inc_atomic": False, "omp": False}
    
Target.register(Seq)
Target.register(Mpi)
Target.register(Omp)
Target.register(Cuda)
Target.register(Hip)
Target.register(Sycl)