#!/usr/bin/env python3

import argparse
import logging
import sys

from sage.cli.eval_cmd import EvalCmd
from sage.cli.interactive_shell_cmd import InteractiveShellCmd
from sage.cli.notebook_cmd import JupyterNotebookCmd
from sage.cli.options import CliOptions
from sage.cli.version_cmd import VersionCmd
from sage.cli.run_file_cmd import RunFileCmd


# Sage CLI options that take a value (i.e. consume the following token, unless
# the value is glued in the same token via "--name=value" or "-xvalue").
_OPTS_WITH_VALUE = ("-c", "--command", "-n", "--notebook")

# Choices for "-n / --notebook".  We need to know these so that we don't
# greedily consume a script filename as the option's value (the option uses
# ``nargs='?'`` and is restricted to these choices).
_NOTEBOOK_CHOICES = ("jupyter", "jupyterlab")


def _split_input_args(input_args):
    """
    Split the raw command-line tokens into

    * ``sage_args``  -- arguments to be parsed by Sage's own ``argparse`` parser
      (this includes the script filename, when there is one, so that it gets
      assigned to the ``file`` positional), and
    * ``script_args`` -- arguments to be forwarded verbatim to the executed
      script (everything after the script filename, or after a literal ``--``).

    The split happens at the first non-option token, which is interpreted as
    the script filename.  Tokens that look like Sage's own options consume
    their value(s) as appropriate so that, e.g., the ``foo.sage`` in
    ``sage -c "1+1" foo.sage`` is recognized as the script.

    This is required because ``argparse`` does not support a "stop here and
    forward the rest" semantics for positional arguments: without this
    pre-split, options such as ``-y`` that the user script wants to handle
    itself would be consumed (and rejected) by Sage's own parser
    (``sage`` issues #40871, #41908).
    """
    sage_args = []
    i = 0
    n = len(input_args)
    while i < n:
        tok = input_args[i]

        # Explicit terminator: everything after `--` is positional.  We hand
        # `--` plus the script filename (if any) to argparse, and forward the
        # remainder to the script.
        if tok == "--":
            sage_args.append(tok)
            i += 1
            if i < n:
                sage_args.append(input_args[i])
                i += 1
            return sage_args, list(input_args[i:])

        # First non-option token: this is the script filename.  Stop here and
        # forward the rest to the script.
        if not tok.startswith("-") or tok == "-":
            sage_args.append(tok)
            return sage_args, list(input_args[i + 1:])

        # Glued long form "--name=value": single token, nothing extra to do.
        if tok.startswith("--") and "=" in tok:
            sage_args.append(tok)
            i += 1
            continue

        # Options that take a value.  Handle short-glued forms like "-cFOO"
        # and the standalone forms "-c FOO" / "--command FOO" / "--command=FOO".
        consumed_value = False
        for opt in _OPTS_WITH_VALUE:
            if tok == opt:
                sage_args.append(tok)
                i += 1
                # `--notebook` / `-n` use ``nargs='?'`` with a fixed set of
                # choices: only consume the next token as the value if it is
                # actually one of those choices, otherwise it is the script.
                if opt in ("-n", "--notebook"):
                    if i < n and input_args[i] in _NOTEBOOK_CHOICES:
                        sage_args.append(input_args[i])
                        i += 1
                else:
                    # `-c` / `--command` use ``nargs='?'``: consume the next
                    # token as the value if it does not look like another
                    # option.  This means `sage -c CMD script.sage` runs CMD
                    # and treats script.sage as the file positional, matching
                    # the prior argparse behavior.
                    if i < n and not input_args[i].startswith("-"):
                        sage_args.append(input_args[i])
                        i += 1
                consumed_value = True
                break
            # short-glued form, e.g. "-cFOO" (only for short options)
            if len(opt) == 2 and tok.startswith(opt) and len(tok) > 2:
                sage_args.append(tok)
                i += 1
                consumed_value = True
                break
        if consumed_value:
            continue

        # Any other dash-prefixed token: a Sage boolean flag, an unknown
        # option, or a combined short flag like ``-vq``.  Hand it to argparse
        # which will either accept it or report a proper error.
        sage_args.append(tok)
        i += 1

    return sage_args, []


def main() -> int:
    input_args = sys.argv[1:]
    parser = argparse.ArgumentParser(
        prog="sage",
        description="If no command is given, starts the interactive interpreter where you can enter statements and expressions, immediately execute them and see their results.",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        default=False,
        help="print additional information",
    )
    parser.add_argument(
        "-q",
        "--quiet",
        action="store_true",
        default=False,
        help="do not display the banner",
    )
    parser.add_argument(
        "--simple-prompt",
        action="store_true",
        default=False,
        help="use simple prompt IPython mode",
    )

    VersionCmd.extend_parser(parser)
    JupyterNotebookCmd.extend_parser(parser)
    EvalCmd.extend_parser(parser)
    RunFileCmd.extend_parser(parser)

    if not input_args:
        return InteractiveShellCmd(CliOptions()).run()

    # Split off arguments destined for a user script *before* invoking the
    # argparse parser, so that flags like `-y` understood by the script are
    # not consumed (or rejected) by Sage's own parser.  See issue #41908.
    sage_args, script_args = _split_input_args(input_args)

    args = parser.parse_args(sage_args)
    options = CliOptions(**vars(args), script_args=script_args)

    logging.basicConfig(level=logging.DEBUG if options.verbose else logging.INFO)

    if args.file:
        return RunFileCmd(options).run()
    if args.command:
        return EvalCmd(options).run()
    if args.notebook:
        return JupyterNotebookCmd(options).run()
    return InteractiveShellCmd(options).run()
