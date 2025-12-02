Original implementations of Superdroplet models in various languages.

These are the implementations I started with (specifically `sd_cython`), which
were forward-ported without necessarily trying to preserve the design across
languages. In November, 2025 I wrote a new, idomatic model from scratch in Rust
and used that to back-port new Python, Modern Fortran, Numba, and Julia implementations
to get better apples-to-apples comparisons.