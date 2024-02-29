#! /bin/sh

julia --project=. benchmark.jl | tee benchmark.md
