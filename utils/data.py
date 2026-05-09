from dataclasses import dataclass, field
from typing import Any


@dataclass
class ResultData:
    label: str
    value: Any
    show: bool = True
    fmt: str = ""

    def __str__(self) -> str:
        v = format(self.value, self.fmt) if self.fmt else str(self.value)
        return f"{self.label}: {v}"

    def __repr__(self) -> str:
        return str(self)
