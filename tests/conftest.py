import sys
from pathlib import Path

# Allow running the test suite without `pip install -e .` first.
SRC = Path(__file__).resolve().parent.parent / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))
