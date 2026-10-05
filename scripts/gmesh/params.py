"""Mesh parameters from <case>/system/gmshParams.

Every key in the file must be known here: a misspelt key is an error, not
a silent fallback to the default.
"""

from __future__ import annotations

import re
from pathlib import Path

DEFAULTS: dict[str, float] = {
    # domain (mm)
    "xInlet": -800.0, "xOutlet": 5600.0, "rFar": 600.0,
    "xZoneStart": -40.0, "xNoseEnd": 40.0, "xAftStart": 620.0, "xZoneEnd": 1200.0,
    "rZone": 120.0,            # circumradius of the polygon P bounding the inner zones
    "coreTip": 0.5, "xTipInterface": -30.0, "coreWake": 20.0, "coreFar": 60.0,
    "xWakeSplit": 820.0,
    # fins
    "xSlabEnd": 900.0,         # end of the fin region (start: xAftStart)
    "finWrap": 6.0,            # hood thickness around each fin
    # spacing (mm)
    "nTheta": 16, "firstLayer": 0.003, "layerRatio": 1.2, "nLayers": 20, "growth": 1.1,
    "firstLayerTip": 0.02,
    "hUpstream": 2.0, "hOgive": 1.0, "hShoulder": 0.5, "hWall": 2.0, "hRingOut": 4.0,
    "hRelax": 2.0, "hWake": 4.0, "hWakeRadial": 1.0, "hSeam": 4.0, "hFarMid": 8.0, "hFar": 30.0,
    "hFinChord": 2.0, "hFinSpan": 1.5, "hSlabOuter": 3.0, "hLE": 0.2,
}
INTEGER = ("nTheta", "nLayers")
_ENTRY = re.compile(r"^\s*([A-Za-z]\w*)\s+([-+0-9.eE]+)\s*;", re.MULTILINE)


def read(case: Path) -> dict[str, float]:
    text = (case / "system" / "gmshParams").read_text()
    text = re.sub(r"FoamFile\s*\{[^}]*\}", "", text)            # header: version 2.0 etc.
    text = re.sub(r"//[^\n]*", "", text)
    found = {k: float(v) for k, v in _ENTRY.findall(text)}
    unknown = sorted(set(found) - set(DEFAULTS))
    if unknown:
        raise ValueError(f"unknown gmshParams keys: {', '.join(unknown)}")
    params = DEFAULTS | found
    for name in INTEGER:
        params[name] = int(params[name])
    return params
