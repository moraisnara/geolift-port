"""
Back-compat shim. The implementation now lives in the installable package
`geolift_fast` (geolift_py/geolift_fast/). This module re-exports it so older
imports (`from market_selection import ...`) keep working. Prefer:

    from geolift_fast.market_selection import Panel, simulate_combo
    # or the high-level API:
    from geolift_fast import power_curves, best_markets
"""
from geolift_fast.market_selection import *          # noqa: F401,F403
from geolift_fast.market_selection import (          # also re-export privates used by scripts
    _stat_func, _scaled_l2,                          # noqa: F401
)
