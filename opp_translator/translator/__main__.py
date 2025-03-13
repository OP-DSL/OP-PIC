import dataclasses
import os
import json
import re
from argparse import ArgumentParser, ArgumentTypeError, Namespace
from datetime import datetime
from pathlib import Path

#custom implementation imports
import cpp
from jinja_utils import env
from language import Lang
from opp import OppError, Type
from scheme import Scheme
from store import Application, ParseError
from target import Target
from util import getVersion, safeFind

def main(argv=None) -> None:

    parser = ArgumentParser(prog="opp-translator")

    parser.add_argument("-V", "--version", help="Version", action="version", version=getVersion()) #this needs version tag
    parser.add_argument("-v", "--verbose", help="Verbose", action="store_false")
    parser.add_argument("-d", "--dump", help="JSON store dump", action="store_true")
    parser.add_argument("-o", "--out", help="Output directory", type=isDirPath)
    parser.add_argument("-c", "--config", help="Target configuration", type=json.loads, default="{}")
    parser.add_argument("-soa", "--force_soa", help="Force structs of arrays", action="store_true")
    parser.add_argument("--suffix", help="Add a suffix to genreated program translations", default="_opp")
    parser.add_argument("-I", help="Add to include directories", type=isDirPath, action="append", nargs=1, default=[])
    parser.add_argument("-i", help="Add to include files", type=isFilePath, action="append", nargs=1, default=[])
    parser.add_argument("-D", help="Add to preprocessor defines", action="append", nargs=1, default=[])
    parser.add_argument("--file_paths", help="Input OPP sources", type=isFilePath, nargs="+")

    target_names = [target.name for target in Target.all()] #TODO: implement Target Findable class
    parser.add_argument("-t", "--target", help="Code-gereration target", type=str, action="append", nargs=1, choices=target_names, default=[])

    #invoking arg parser
    args = parser.parse_args(argv)

    if os.environ.get("OPP_AUTO_SOA") is not None:
        args.force_soa = True

    file_parents = [Path(file_path).parent for file_path in args.file_paths]

    if args.out is None:
        args.out = file_parents[0]

    # Generating code for targets
    for [target] in args.target:
        print(f'1 target {Target.find(target)}')

    #checking includes of OPP
    if os.environ.get("OPP_INSTALL_PATH") is not None:
        ops_install_path = Path(os.environ.get("OPP_INSTALL_PATH"))
        args.I = [[str(ops_install_path/"include")]] + args.I
    else:
        script_parents = list(Path(__file__).resolve().parents)
        if len(script_parents) >= 3 and script_parents[2].stem == "OPP":
            args.I = [[str(script_parents[2].joinpath("opp_lib/include"))]] + args.I

    for [target] in args.target:
        print(f'2 target {Target.find(target)}')

    args.I = [[str(file_parent)] for file_parent in file_parents] + args.I

    # Collect the set of file extensions 
    extensions = {str(Path(file_path).suffix)[1:] for file_path in args.file_paths}

    if not extensions:
        exit("Missing file extensions, unable to determine target language.")
    elif len(extensions) > 1:
        exit("Varying file extensions, unable to determine target language.")
    else:
        [extension] = extensions

    lang = Lang.find(extension)

    if lang is None:
        exit(f"Unknown file extension: {extension}")

    Type.set_formatter(lang.formatType)

    if len(args.target) == 0:
        args.target = [[target_name] for target_name in target_names]

    print("Code-gen : Program parsing phase started......")
    try:
        app = parse(args, lang)
    except ParseError as e:
        exit(e)

    if args.verbose:
        print()
        print(app)

    try: 
        validate(args, lang, app)
    except OppError as e:
        exit(e)

    print("Code-gen : Program parsing phase finished.........")

    # Generate program translations
    print("Code-gen : Program translation phase started......")
    app_consts = app.consts()
    for i, program in enumerate(app.programs, 1):
        include_dirs = set([Path(dir) for [dir] in args.I])
        defines = [define for [define] in args.D]

        source = lang.translateProgram(program, include_dirs, defines, app_consts, args.force_soa)
        
        if not args.force_soa and program.soa_val:
            args.force_soa = program.soa_val

        new_file = os.path.splitext(os.path.basename(program.path))[0]
        ext = os.path.splitext(os.path.basename(program.path))[1]
        new_path = Path(args.out, f"{new_file}{args.suffix}{ext}")

        with open(new_path, "w") as new_file:
            new_file.write(f"\n{lang.com_delim} Auto-generated at {datetime.now()} by opp-translator\n")
            new_file.write(source)

            if args.verbose:
                print(f"Translated program {i} of {len(args.file_paths)}: {new_path}")
    
    print("Code-gen : Program translation phase finished.........")

    # Generating code for targets
    for [target] in args.target:
        target = Target.find(target)

        # Applying user defined configs to the target config
        for key in target.config:
            if key in args.config:
                target.config[key] = args.config[key]

        scheme = Scheme.find((lang, target))

        if not scheme:
            if args.verbose:
                print(f"No scheme register for {lang}/{target}")

            continue

        if args.verbose:
            print(f"Translation scheme: {scheme}")

        print("Code-gen : Generating target specific template, scheme - " + scheme.target.name)
        codegen(args, scheme, app, args.force_soa)

        if args.verbose:
            print(f"Translation completed: {scheme}")


def parse(args: Namespace, lang: Lang) -> Application:
    app = Application()

    # Collect the include directories
    include_dirs = set([Path(dir) for [dir] in args.I])
    defines = [define for [define] in args.D]

    # Parse the input files
    for i, raw_path in enumerate(args.file_paths, 1):
        if args.verbose:
            print(f"Parsing file {i} of {len(args.file_paths)}: {raw_path}")

        # Parse the program
        program = lang.parseProgram(Path(raw_path), include_dirs, defines)
        app.programs.append(program)

    return app

def validate(args: Namespace, lang: Lang, app: Application) -> None:
    # Run sementic checks on the application
    app.validate(Lang)

    if args.dump:
        store_path = Path(args.out, "store.json")
        serializer = lambda o: getattr(o, "__dict__", "unserializable")

        # Write application dump
        with open(store_path, "w") as file:
            file.write(json.dumps(app, default=serializer, indent=4))

        if args.verbose:
            print("Dumped store: ", store_path.resolve(), end="\n\n")

def write_file(path: Path, text: str) -> None:
    if path.is_file():
        prev_text = path.read_text()

        if text == prev_text:
            return

    with path.open("w") as f:
        # f.write(f"{scheme.lang.com_delim} Auto-generated at {datetime.now()} by op2-translator\n\n")
        f.write(text)

def codegen(args: Namespace, scheme: Scheme, app: Application, force_soa: bool = False) -> None:
    # Collect the paths of the generated files
    include_dirs = set([Path(dir) for [dir] in args.I])
    defines = [define for [define] in args.D]

    # Generate loop hosts
    for i, (loop, program) in enumerate(app.uniqueLoops(), 1):
        
        force_generate = scheme.target == Target.find("seq")

        # Generate loop host source
        res = scheme.genLoopHost(env, loop, program, app, i, args.config, force_generate)

        if res is None:
            print(f"Error: unable to generate loop host {i}")
            continue

        source, extension, fallback = res

        # Form output file path
        Path(args.out, scheme.target.name).mkdir(parents=True, exist_ok=True)
        path = Path(
            args.out,
            scheme.target.name,
            f"{loop.kernel}_loop.{extension}",
        )

        # Write the generated source file
        write_file(path, source)

        if args.verbose and not fallback:
            print(f"Generated loop host {i} of {len(app.loops())}: {path}")

        if args.verbose and fallback:
            print(f"Generated loop host {i} of {len(app.loops())} (fallback): {path}")


    # Generate master kernel file
    if scheme.master_kernel_template is not None:
        user_types_name = f"user_types.{scheme.lang.include_ext}"
        user_types_candidates = [Path(dir, user_types_name) for dir in include_dirs]
        user_types_file = safeFind(user_types_candidates, lambda p: p.is_file())

        source, name = scheme.genMasterKernel(env, app, user_types_file)

        Path(args.out, scheme.target.name).mkdir(parents=True, exist_ok=True)
        path = Path(args.out, scheme.target.name, name)

        write_file(path, source)

        if args.verbose:
            print(f"Generated master kernel file: {path}")


def isDirPath(path):
    if os.path.isdir(path):
        return path
    else:
        raise ArgumentTypeError("Invalid directory path: {path}")

def isFilePath(path):
    if os.path.isfile(path):
        return path
    else:
        raise ArgumentTypeError("Invalid file: {path}")

def print_ast(cursor, depth=0):
    indent = '  ' * depth
    print(f"{indent}{cursor.kind.name}: {cursor.spelling}")
    for child in cursor.get_children():
        print_ast(child, depth + 1)
    

if __name__ == "__main__":
    main()
