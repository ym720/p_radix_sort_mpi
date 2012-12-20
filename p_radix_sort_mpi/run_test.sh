#!/bin/bash

SIZES="1000 10000 100000 1000000 10000000"

for SIZE in $SIZES; do
  echo "=== radix_sort ==="
  ./radix_sort numbers.in $SIZE
  sleep 1
  echo ""
  echo "=== p_radix_sort ==="
  mpiexec -n 2 p_radix_sort numbers.in $SIZE
  sleep 1
  echo ""
done

