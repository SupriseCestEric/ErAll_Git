#!/usr/bin/bash

find . -name '*:*' -type f -print0 | perl -0ne ' rename $_, s{[^/]+$}{$& =~ y/:/-/r}res or warn "rename $_: $!"'
