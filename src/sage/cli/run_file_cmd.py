import argparse
import sys
import types

from sage.cli.options import CliOptions
from sage.repl.preparse import preparse_file_named
from sage.repl.load import load_cython
from sage.misc.temporary_file import tmp_filename
from sage.all import sage_globals


class RunFileCmd:
    @staticmethod
    def extend_parser(parser: argparse.ArgumentParser):
        r"""
        Extend the parser with the "run file" command.

        INPUT:

        - ``parsers`` -- the parsers to extend.

        OUTPUT:

        - the extended parser.
        """
        parser.add_argument(
            "file",
            nargs="*",
            help="execute the given file as sage code",
        )

    def __init__(self, options: CliOptions):
        r"""
        Initialize the command.
        """
        self.options = options
        # Rebuild ``sys.argv`` so that the executed script sees the same
        # arguments it would receive under plain Python:
        # ``sys.argv == [<script>, *script_args]``.  ``script_args`` is the
        # tail of the original command line that was split off from Sage's
        # own arguments by ``sage.cli._split_input_args`` (see issue #41908).
        script = options.file[0] if options.file else sys.argv[0]
        forwarded = list(options.script_args or [])
        sys.argv = [script, *forwarded]

    def run(self) -> int:
        r"""
        Execute the given command.
        """
        source_file = self.options.file[0]
        input_file = source_file
        if input_file.endswith('.sage'):
            input_file = str(preparse_file_named(input_file))

        main_module = types.ModuleType('__main__')
        main_module.__dict__.update(sage_globals())
        main_module.__name__ = '__main__'
        main_module.__file__ = source_file
        main_module.__package__ = None
        main_module.__loader__ = None
        main_module.__spec__ = None
        main_module.__cached__ = None
        main_module.__doc__ = None

        previous_main = sys.modules.get('__main__')
        sys.modules['__main__'] = main_module
        try:
            if input_file.endswith('.pyx') or input_file.endswith('.spyx'):
                s = load_cython(input_file)
                exec(compile(s, tmp_filename(), 'exec'), main_module.__dict__)
            else:
                with open(input_file, 'rb') as f:
                    source = f.read()
                exec(compile(source, input_file, 'exec'), main_module.__dict__)
        finally:
            if previous_main is None:
                del sys.modules['__main__']
            else:
                sys.modules['__main__'] = previous_main
        return 0
