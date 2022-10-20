import sys
import pickle
import base64
import traceback
import importlib.util
from typing import Callable, Any


def wildcardImportFromFile(path: str) -> None:
    """ same as from main.py import * """
    spec = importlib.util.spec_from_file_location("tmp", path)
    assert spec
    assert spec.loader
    module = importlib.util.module_from_spec(spec)
    assert module
    spec.loader.exec_module(module)
    for i in dir(module):
        if not i.startswith("__"):
            globals()[i] = getattr(module, i)


def isolateRun(func: Callable, input_name: str) -> Any:
    """ run the function """
    try:
        return func(input_name)
    except BaseException as e:
        print(e)
        return BaseException(traceback.format_exc())


if __name__ == "__main__":
    wildcardImportFromFile(sys.argv[1])
    name = sys.argv[2]
    func, input_name = pickle.load(open(f"{name}.tmp.in", "rb"))
    out = isolateRun(func, input_name)
    pickle.dump(out, open(f"{name}.tmp.out", "wb"), pickle.HIGHEST_PROTOCOL)
