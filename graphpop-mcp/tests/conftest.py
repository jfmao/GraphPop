"""Add project root to path so tooluniverse/ is importable in tests."""

import sys
from pathlib import Path

# Add the graphpop-mcp root so 'tooluniverse' package is importable
project_root = Path(__file__).parent.parent
if str(project_root) not in sys.path:
    sys.path.insert(0, str(project_root))
