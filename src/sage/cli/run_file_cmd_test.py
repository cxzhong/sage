import argparse
import sys
from unittest.mock import patch

from sage.cli.run_file_cmd import RunFileCmd
from sage.cli.options import CliOptions


def _parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("-q", "--quiet", action="store_true", default=False)
    RunFileCmd.extend_parser(parser)
    return parser


def test_run_file_parser_forwards_arguments_after_file():
    args = _parser().parse_args(["-q", "test.sage", "42", "-y", "7", "--quiet"])

    assert args.quiet is True
    assert args.file == ["test.sage", "42", "-y", "7", "--quiet"]


def test_run_file_parser_handles_explicit_terminator():
    args = _parser().parse_args(["--", "-test.sage", "-y"])

    assert args.file == ["-test.sage", "-y"]


def test_run_file_cmd(capsys, tmp_path):
    file = tmp_path / "test.sage"
    file.write_text("print(3^33)")
    options = CliOptions(file=[str(file)])
    run_file_cmd = RunFileCmd(options)

    run_file_cmd.run()
    captured = capsys.readouterr()
    assert captured.out == "5559060566555523\n"


def test_run_file_cmd_with_args(capsys, tmp_path):
    with patch.object(sys, 'argv', ["python3", "test.sage", "1", "1"]):
        file = tmp_path / "test.sage"
        file.write_text("import sys; print(int(sys.argv[1]) + int(sys.argv[2]))")
        options = CliOptions(file=[str(file), "1", "1"])
        run_file_cmd = RunFileCmd(options)

        run_file_cmd.run()
        captured = capsys.readouterr()
        assert captured.out == "2\n"


def test_run_file_cmd_argv_matches_python(capsys, tmp_path):
    """``sys.argv`` inside a script must look like ``[script, *args]``."""
    with patch.object(sys, 'argv', ["python3", "test.sage", "-y", "7"]):
        file = tmp_path / "test.sage"
        file.write_text("import sys; print(sys.argv)")
        options = CliOptions(file=[str(file), "-y", "7"])
        RunFileCmd(options).run()
        captured = capsys.readouterr()
        assert captured.out == f"[{str(file)!r}, '-y', '7']\n"


def test_run_file_cmd_dunder_name_is_main(capsys, tmp_path):
    """Inside the script, ``__name__`` must equal ``'__main__'`` (#42159)."""
    file = tmp_path / "test.sage"
    file.write_text(
        "import __main__\n"
        "print(__name__)\n"
        "print(__file__)\n"
        "print(__package__)\n"
        "print(__spec__)\n"
        "print(__main__.__name__)\n"
        "print(__main__.__file__)\n"
    )
    options = CliOptions(file=[str(file)])
    RunFileCmd(options).run()
    captured = capsys.readouterr()
    assert captured.out == (
        "__main__\n"
        f"{file}\n"
        "None\n"
        "None\n"
        "__main__\n"
        f"{file}\n"
    )


def test_run_file_cmd_reset_preserves_dunder_name(capsys, tmp_path):
    """``reset()`` must not clobber the script's ``__main__`` metadata."""
    file = tmp_path / "test.sage"
    file.write_text(
        "import __main__\n"
        "print(__name__)\n"
        "print(__file__)\n"
        "print(__package__)\n"
        "print(__spec__)\n"
        "reset()\n"
        "print(__name__)\n"
        "print(__file__)\n"
        "print(__package__)\n"
        "print(__spec__)\n"
        "print(__main__.__name__)\n"
        "print(__main__.__file__)\n"
    )
    options = CliOptions(file=[str(file)])
    RunFileCmd(options).run()
    captured = capsys.readouterr()
    assert captured.out == (
        "__main__\n"
        f"{file}\n"
        "None\n"
        "None\n"
        "__main__\n"
        f"{file}\n"
        "None\n"
        "None\n"
        "__main__\n"
        f"{file}\n"
    )
