"""
Time tracking utilities.
"""

import dataclasses
import functools


@dataclasses.dataclass
@functools.total_ordering
class Stopwatch:
    """A stopwatch to measure elapsed time in minutes and seconds.

    Attributes:
        minutes: The number of minutes elapsed.
        seconds: The number of seconds elapsed.

    Example:
        >>> stopwatch = Stopwatch(minutes=1, seconds=30)
        >>> stopwatch.increment(10)
        >>> print(stopwatch)
        Stopwatch(minutes=1, seconds=40)
    """

    minutes: int = dataclasses.field(default=0)
    seconds: int = dataclasses.field(default=0)

    def __init__(self, seconds: int = 0):
        seconds_over = seconds % 60
        self.minutes = (seconds - seconds_over) // 60
        self.seconds = seconds_over

    def increment(self, seconds: int) -> None:
        """Increment the stopwatch by the given number of seconds."""
        self.seconds += seconds
        if self.seconds >= 60:
            seconds_overflow = self.seconds % 60
            minutes_overflow = (self.seconds - seconds_overflow) // 60
            self.minutes += minutes_overflow
            self.seconds = seconds_overflow

    def total_seconds(self) -> int:
        """Return the total number of seconds elapsed."""
        return self.minutes * 60 + self.seconds

    def __str__(self) -> str:
        return f"{self.minutes:4d} min {self.seconds:2d} sec"

    def __repr__(self) -> str:
        return f"Stopwatch(minutes={self.minutes}, seconds={self.seconds})"

    def __eq__(self, other: "Stopwatch") -> bool:
        return self.minutes == other.minutes and self.seconds == other.seconds

    def __lt__(self, other: "Stopwatch") -> bool:
        return self.total_seconds() < other.total_seconds()
