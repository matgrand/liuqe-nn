"""Automatically discover and expose dissertation figure modules."""
from importlib import import_module
from pkgutil import iter_modules

from . import _common
from ._common import compile_inspection_pdf


def configure(SHOW_PLOTS=True, WAIT_FOR_PLOTS=False):
    """Configure shared display behavior for all figures."""
    _common.SHOW_PLOTS = SHOW_PLOTS
    _common.WAIT_FOR_PLOTS = WAIT_FOR_PLOTS


__all__ = ['compile_inspection_pdf', 'configure']
for _module_info in iter_modules(__path__):
    if _module_info.name.startswith('_') or _module_info.name in ('explore_raw_ds', 'find_shapes'): continue
    globals()[_module_info.name] = import_module(f'{__name__}.{_module_info.name}')
    __all__.append(_module_info.name)
