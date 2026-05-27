# Python BEM - Blade Element Momentum Theory Software.

# Copyright (C) 2022 M. Smrekar

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.append(str(ROOT))

from calculation_runner import calculate_power_3d


def test_calculate_power_3d_from_sample_file():
    settings_path = ROOT / "karlsen.bem"
    with settings_path.open() as f:
        settings = json.load(f)

    results = calculate_power_3d(settings, print_progress=False)

    assert isinstance(results, dict)
    assert "power" in results
    assert len(results["power"]) > 0
