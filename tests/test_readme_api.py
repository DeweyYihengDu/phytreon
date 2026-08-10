"""The README's code examples have to refer to APIs that actually exist.

The README is the entry point most people read first, and nothing else checks
it -- a renamed parameter or a dropped function stays documented indefinitely,
and the reader finds out by getting a TypeError. This walks every ``pt.*`` call
in every fenced python block and verifies the attribute chain resolves and each
keyword argument is a real parameter.

It cannot check that the examples *run* -- they use placeholder names like
``asvs.fasta`` and ``trait`` on purpose -- but "this function takes that
argument" is exactly the part that goes stale silently.
"""
import ast
import inspect
import pathlib
import re

import matplotlib
matplotlib.use("Agg")

import phytreon as pt

README = pathlib.Path(__file__).resolve().parent.parent / "README.md"


def _python_blocks():
    md = README.read_text(encoding="utf-8")
    return re.findall(r"```python\n(.*?)```", md, re.S)


def test_readme_has_python_examples_to_check():
    # if the extraction regex ever stops matching, every other test in this file
    # would pass by vacuously checking nothing
    blocks = _python_blocks()
    assert len(blocks) >= 10
    assert sum("pt." in b for b in blocks) >= 8


def test_every_readme_python_block_parses():
    # a block that does not parse is a block nobody can copy and paste, and it is
    # also a block the API check below skips in silence -- which is how two
    # unparseable blocks (a bare `...` inside a dict literal) once hid the entire
    # set of examples for a release's new functions
    problems = []
    for i, block in enumerate(_python_blocks()):
        try:
            ast.parse(block)
        except SyntaxError as exc:
            problems.append(f"block {i}: {exc}")
    assert not problems, "unparseable README example(s): " + "; ".join(problems)


def _pt_calls():
    """Every ``pt.<...>(...)`` call in the README, as (block index, path, node)."""
    for bi, block in enumerate(_python_blocks()):
        try:
            tree = ast.parse(block)
        except SyntaxError:
            continue          # reported by the test above; not this one's job
        for node in ast.walk(tree):
            if not isinstance(node, ast.Call):
                continue
            parts, cur = [], node.func
            while isinstance(cur, ast.Attribute):
                parts.append(cur.attr)
                cur = cur.value
            if isinstance(cur, ast.Name) and cur.id == "pt":
                yield bi, list(reversed(parts)), node


def test_every_documented_api_exists():
    missing = []
    for bi, parts, _ in _pt_calls():
        obj, path = pt, "pt"
        for p in parts:
            path += "." + p
            if not hasattr(obj, p):
                missing.append(f"block {bi}: {path}")
                break
            obj = getattr(obj, p)
    assert not missing, "README references APIs that do not exist: " + "; ".join(missing)


def test_every_documented_keyword_argument_is_a_real_parameter():
    wrong = []
    checked = 0
    for bi, parts, node in _pt_calls():
        obj = pt
        for p in parts:
            if not hasattr(obj, p):
                obj = None
                break
            obj = getattr(obj, p)
        if obj is None or not callable(obj):
            continue
        try:
            sig = inspect.signature(obj)
        except (ValueError, TypeError):
            continue          # builtins and C-level callables have no signature
        checked += 1
        names = set(sig.parameters)
        takes_kwargs = any(p.kind is p.VAR_KEYWORD
                           for p in sig.parameters.values())
        for kw in node.keywords:
            if kw.arg and kw.arg not in names and not takes_kwargs:
                wrong.append(
                    f"block {bi}: pt.{'.'.join(parts)}() has no parameter "
                    f"{kw.arg!r}"
                )
    assert checked >= 20, f"only {checked} calls resolved -- extraction may be broken"
    assert not wrong, "README passes arguments that do not exist: " + "; ".join(wrong)
