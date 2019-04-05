# nim-falcon
Nim versions of FALCON executables

See https://github.com/PacificBiosciences/FALCON

Requires nim >= 0.19.9

## Installation
Our packages are in flux, not officially "Nimble" yet.

So do this:
```
export NIMBLE_DIR=~/.nimble
git submodule update --init --recursive
make install-all

nimble install

./link-exe.sh ${PREFIX}
```

The executables lack the `.exe` that our workflow expects, so
we use `link-exe.sh` to rename them. But that's the easy part.
If you get everything built, you're almost home.

## Goals (eventually)
* A single executable, or very small number of them
* On Linux, a static binary
