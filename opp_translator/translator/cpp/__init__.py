import os
from functools import lru_cache
from io import StringIO
from pathlib import Path
from typing import FrozenSet, List, Optional, Set, Tuple, Any

import clang.cindex

import pcpp

import cpp.parser
import cpp.translator.program
import opp as OPP
from language import Lang
from store import Application, Location, ParseError, Program

libclang_path = os.getenv("LIBCLANG_PATH")
if libclang_path is not None:
    clang.cindex.Config.set_library_file(libclang_path)

class Preprocessor(pcpp.Preprocessor):
    def __init__(self, lexer=None):
        super(Preprocessor, self).__init__(lexer)
        self.line_directive = None

    def on_comment(self, tok: str) -> bool:
        return True

    def on_error(self, file: str, line: int, msg: str) -> None:
        loc = Location(file, line, 0)
        raise ParseError(msg, loc)

    def on_include_not_found(self, is_malformed, is_system_include, curdir, includepath) -> None:
        if is_system_include:
            raise pcpp.OutputDirective(pcpp.Action.IgnoreAndPassThrough)

        super().on_include_not_found(is_malformed, is_system_include, curdir, includepath)

def print_ast(cursor, depth=0):
    indent = '  ' * depth
    print(f"{indent}{cursor.kind.name}: {cursor.spelling}")
    for child in cursor.get_children():
        print_ast(child, depth + 1)

class Cpp(Lang):
    name = "C++"

    source_exts = ["cpp"]
    include_ext = "h"
    kernel_dir = True

    com_delim = "//"
    zero_idx = True

    @lru_cache(maxsize=None)
    def parseFile(
        self, path: Path, include_dirs: FrozenSet[Path], defines: FrozenSet[str], preprocess: bool = False
        ) -> Tuple[clang.cindex.TranslationUnit, str]:
        args = [f"-I{dir}" for dir in include_dirs]
        args = args + [f"-D{define}" for define in defines]
        args = args +['-std=c++11']
        source = path.read_text()

        if preprocess:
            preprocessor = Preprocessor() 

            for dir in include_dirs:
                preprocessor.add_path(str(dir.resolve()))

            for define in defines:
                if "=" not in define:
                    define = f"{define}=1"

                preprocessor.define(define.replace("=", " ", 1))

            preprocessor.parse(source, str(path.resolve()))
            source_io = StringIO()
            preprocessor.write(source_io)

            source_io.seek(0)
            source = source_io.read()

        translation_unit = clang.cindex.Index.create().parse(
            path,
            unsaved_files=[(path, source)],
            args=args,
            options=clang.cindex.TranslationUnit.PARSE_DETAILED_PROCESSING_RECORD
        )

        for diagnostic in iter(translation_unit.diagnostics):
            if not "'stddef.h' file not found" in diagnostic.spelling:
                print(
                    f"Clang parse diagnotic message: severity {diagnostic.severity} at: "
                    f"{cpp.parser.parseLocation(diagnostic)}: {diagnostic.spelling}"
                )

        return translation_unit, source

    def parseProgram(self, path: Path, include_dirs: Set[Path], defines: List[str]) -> Program:
        print("Code-gen, parsing file:  " + str(path))
        ast, source = self.parseFile(path, frozenset(include_dirs), frozenset(defines))
        ast_pp, source_pp =  self.parseFile(path, frozenset(include_dirs), frozenset(defines), preprocess = True)

        # print_ast(ast.cursor)

        program = Program(path, ast_pp, source_pp)

        cpp.parser.parseLoops(ast, program)
        cpp.parser.parseMeta(ast_pp.cursor, program)

        program.enrichLoopData()

        return program

    def translateProgram(self, program: Program, include_dirs: Set[Path], defines: List[str], app_consts: List[OPP.Const], force_soa: bool = False) -> str:
        return cpp.translator.program.translateProgram(program.path.read_text(), program, force_soa)

    def formatType(self, typ: OPP.Type) -> str:
        int_types = {
            (True, 16): "short",
            (True, 32): "int",
            (True, 64): "long long",
            (False, 16): "unsigned short",
            (False, 32): "unsigned int",
            (False, 64): "unsigned long long",
        }

        float_type = {16: "half", 32: "float", 64: "double"}

        if isinstance(typ, OPP.Int):
            return int_types[(typ.signed, typ.size)]
        elif isinstance(typ, OPP.Float):
            return float_type[typ.size]
        elif isinstance(typ, OPP.Bool):
            return "bool"
        elif isinstance(typ, OPP.Char):
            return "char"
        elif isinstance(typ, OPP.ComplexD):
            return "complexd"
        elif isinstance(typ, OPP.ComplexF):
            return "complexf"
        elif isinstance(typ, OPP.Custom):
            return typ.name
        else:
            assert False

Lang.register(Cpp)

import cpp.schemes
