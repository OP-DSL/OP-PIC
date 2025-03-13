from __future__ import annotations

import sys
import traceback
import re
from abc import abstractmethod
from pathlib import Path 
from typing import Any, Dict, List, Optional, Set, Tuple

from jinja2 import Environment

import opp as OPP
from language import Lang
from store import Application, Program
from target import Target
from util import sycl_set_flat_parallel
from util import extract_intrinsic_functions
from util import Findable
from util import KernelProcess

class Scheme(Findable):
    lang: Lang
    target: Target

    loop_host_template: Path
    move_loop_host_template: Path
    master_kernel_template: Optional[Path]

    def __str__(self) -> str:
        return f"{self.lang.name}/{self.target.name}"

    def canGenLoopHost(self, loop: OPP.Loop) -> bool:
        return True

    def getBaseConfig(self, loop: OPP.Loop) -> Dict[str, Any]:
        return self.target.defaultConfig()

    def getConfig(self, loop: OPP.Loop, config_overrides: List[Dict[str, Dict[str, Any]]]) -> Dict[str, Any]:
        config = self.getBaseConfig(loop)

        for override in config_overrides:
            for loop_match, items in override.items():
                if re.match(loop_match, loop.name):
                    config |= items

        return config
        
    def genLoopHost(
        self,
        env: Environment,
        loop: OPP.Loop,
        program: Program,
        app: Application,
        kernel_idx: int,
        config_overrides: List[Dict[str, Dict[str, Any]]],
        force_generate: bool = False,
    ) -> Optional[Tuple[str, str, bool]]:
        
        if loop.loop_type == OPP.LoopType.PAR_LOOP:
            template = env.get_template(str(self.loop_host_template))
            extension = self.loop_host_template.suffixes[-2][1:]
        else:
            template = env.get_template(str(self.move_loop_host_template))
            extension = self.move_loop_host_template.suffixes[-2][1:]            
        # Load the loop host template

        args = {
            "OP": OPP,
            "lh": loop,
            "kernel_idx": kernel_idx,
            "lang": self.lang,
            "config": self.getConfig(loop, config_overrides),
            "kernel_func": None,
            "host_kernel_func": None,
        }

        try:
            if (not loop.fallback and self.canGenLoopHost(loop)) or force_generate:
                args["host_kernel_func"] = self.translateHostKernel(loop, program, app, args["config"], kernel_idx)
                args["kernel_func"] = self.translateKernel(loop, program, app, args["config"], kernel_idx)

        except Exception as e:
            print(f"Error: kernel translation for kernel {kernel_idx} ({loop.name}) failed ({self}):")
            print(f"  fallback: {loop.fallback}, can generate: {self.canGenLoopHost(loop)}, force_generate: {force_generate}")
            traceback.print_exc(file=sys.stdout)

        if args["kernel_func"] is None and self.fallback is None:
            # print('if args["kernel_func"] is None and self.fallback is None')
            return None

        if self.fallback is not None:
            fallback_wrapper_template = env.get_template(str(self.lang.fallback_wrapper_template))
            fallback_template = env.get_template(str(self.fallback.loop_host_template))

            fallback_args = dict(args)

            fallback_args["config"] = self.fallback.getConfig(loop, config_overrides)
            fallback_args["kernel_func"] = self.fallback.translateKernel(
                loop, program, app, fallback_args["config"], kernel_idx
            )

        if args["kernel_func"] is None:
            return (fallback_template.render(**fallback_args, variant=""), extension, True)

        if self.fallback is None:
            return (template.render(**args, variant=""), extension, False)

        source = template.render(**args, variant="_main")

        source += "\n\n"
        source += fallback_template.render(**fallback_args, variant="_fallback")

        source += "\n\n"
        source += fallback_wrapper_template.render(**args)

        return (source, extension, False)

    def genMasterKernel(self, env: Environment, app: Application, user_types_file: Optional[Path]) -> Tuple[str, str]:
        if self.master_kernel_template is None:
            exit(f"No master kernel template registered for {self}")

        user_types = None
        if user_types_file is not None:
            user_types = user_types_file.read_text()

        # Load the loop host template
        template = env.get_template(str(self.master_kernel_template))

        extension = self.master_kernel_template.suffixes[-2][1:]
        name = f"opp_kernels.{extension}"

        unique_loops = []
        for loop, _ in app.loops():
            if loop.kernel not in unique_loops:
                unique_loops.append(loop.kernel)

        # print(f'at GEN MASTER KERNEL {app.consts()}')
        # Generate source from the template
        return template.render(OP=OPP, app=app, lang=self.lang, target=self.target, user_types=user_types, header=self.header, unique_loops=unique_loops), name

        
    @abstractmethod
    def translateKernel(
        self,
        loop: OPP.Loop,
        program: Program,
        app: Application,
        config: Dict[str, Any],
        kernel_idx: int,
    ) -> str:
        pass
    
    def translateHostKernel(
        self,
        loop: OPP.Loop,
        program: Program,
        app: Application,
        config: Dict[str, Any],
        kernel_idx: int,
    ) -> str:
        return "// Host kernel not used"

    def matches(self, key: Tuple[Lang, Target]) -> bool:
        if not (isinstance(key, tuple) and len(key) == 2 and isinstance(key[0], Lang) and isinstance(key[1], Target)):
            return False

        return self.lang == key[0] and self.target == key[1]