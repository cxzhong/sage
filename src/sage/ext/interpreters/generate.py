"""Generate fast-callable interpreters with compatibility fixes.

Some distribution builds generate these sources with the released version of
``sage_setup``.  Apply the changes needed by this version of Sage when that
older generator is used.  The replacements are no-ops for the in-tree
generator.
"""

import subprocess
import sys
from pathlib import Path


def _insert_declaration(path, marker, declaration):
    text = path.read_text()
    if declaration in text:
        return
    if text.count(marker) != 1:
        raise RuntimeError(f"cannot locate declaration marker in {path}")
    path.write_text(text.replace(marker, marker + "\n" + declaration + "\n", 1))


def _replace_legacy(path, legacy, replacement):
    text = path.read_text()
    if replacement in text:
        return
    if text.count(legacy) != 1:
        raise RuntimeError(f"cannot locate legacy code in {path}")
    path.write_text(text.replace(legacy, replacement, 1))


def _patch_generated_sources(output_dir):
    _insert_declaration(
        output_dir / "interp_cc.c",
        "#include <mpc.h>\n",
        "int cc_py_call_helper(PyObject *, PyObject *, int, mpc_t *, __mpc_struct *);",
    )
    _insert_declaration(
        output_dir / "interp_cdf.c",
        "#include <complex.h>\n",
        "int cdf_py_call_helper(PyObject *, int, double complex *, double complex *);",
    )
    _insert_declaration(
        output_dir / "interp_el.c",
        "#include <Python.h>\n",
        "PyObject *el_check_element(PyObject *, PyObject *);",
    )
    _insert_declaration(
        output_dir / "interp_rr.c",
        "#include <mpfr.h>\n",
        "int rr_py_call_helper(PyObject *, PyObject *, int, mpfr_t *, __mpfr_struct *);",
    )

    interp_cc = output_dir / "interp_cc.c"
    _replace_legacy(
        interp_cc,
        "        mpc_t retval,\n",
        "        mpfr_ptr retval_re, mpfr_ptr retval_im,\n",
    )
    _replace_legacy(
        interp_cc,
        "        mpc_set(retval, i0, MPC_RNDNN);\n",
        "        mpfr_set(retval_re, mpc_realref(i0), MPFR_RNDN);\n"
        "        mpfr_set(retval_im, mpc_imagref(i0), MPFR_RNDN);\n",
    )

    wrapper_cc = output_dir / "wrapper_cc.pyx"
    _replace_legacy(
        wrapper_cc,
        "        mpc_t retval,\n        mpc_t* constants,\n",
        "        mpfr_ptr retval_re, mpfr_ptr retval_im,\n        mpc_t* constants,\n",
    )
    _replace_legacy(
        wrapper_cc,
        "interp_cc(c_args\n            , (<mpc_t>(retval.__re))\n",
        "interp_cc(c_args\n            , retval.__re, retval.__im\n",
    )
    _replace_legacy(
        wrapper_cc,
        "interp_cc(args\n            , result\n",
        "interp_cc(args\n            , mpc_realref(result), mpc_imagref(result)\n",
    )


def main():
    if len(sys.argv) != 3:
        raise SystemExit(f"usage: {sys.argv[0]} GENERATOR OUTPUT_DIR")
    generator = Path(sys.argv[1]).resolve()
    output_dir = Path(sys.argv[2]).resolve()
    # Meson only invokes this command when its declared outputs are stale.
    # Remove the generator's timestamp sentinel so that a partially cleaned
    # output set is regenerated in full.
    (output_dir / "all.py").unlink(missing_ok=True)
    subprocess.run([sys.executable, generator, output_dir], check=True)
    _patch_generated_sources(output_dir)


if __name__ == "__main__":
    main()
