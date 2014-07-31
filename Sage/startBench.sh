#!/bin/sh
sage <<EOF
run "field_tests/mullenTable_benchmark.sage"
%time main()
exit()
EOF
