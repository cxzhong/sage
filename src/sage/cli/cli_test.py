from sage.cli import _split_input_args


def test_split_input_args_forwards_options_after_script_filename():
    sage_args, script_args = _split_input_args(["script.sage", "42", "-y", "7"])

    assert sage_args == ["script.sage"]
    assert script_args == ["42", "-y", "7"]


def test_split_input_args_keeps_sage_options_before_script_filename():
    sage_args, script_args = _split_input_args(["-v", "script.sage", "42", "-q"])

    assert sage_args == ["-v", "script.sage"]
    assert script_args == ["42", "-q"]


def test_split_input_args_handles_explicit_terminator():
    sage_args, script_args = _split_input_args(["--", "-script.sage", "-y"])

    assert sage_args == ["--", "-script.sage"]
    assert script_args == ["-y"]


def test_split_input_args_keeps_command_option_values():
    assert _split_input_args(["-c", "print(2+2)"]) == (["-c", "print(2+2)"], [])
    assert _split_input_args(["--command=print(3+4)"]) == (["--command=print(3+4)"], [])
    assert _split_input_args(["-cprint(5+6)"]) == (["-cprint(5+6)"], [])


def test_split_input_args_keeps_notebook_choice_values():
    assert _split_input_args(["-n", "jupyterlab"]) == (["-n", "jupyterlab"], [])
    assert _split_input_args(["--notebook=jupyterlab"]) == (["--notebook=jupyterlab"], [])


def test_split_input_args_forwards_sage_option_names_after_file():
    sage_args, script_args = _split_input_args(
        ["script.sage", "-c", "1", "-n", "jupyterlab", "-v", "-q", "--simple-prompt"]
    )

    assert sage_args == ["script.sage"]
    assert script_args == ["-c", "1", "-n", "jupyterlab", "-v", "-q", "--simple-prompt"]