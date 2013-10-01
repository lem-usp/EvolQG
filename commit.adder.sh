#!/bin/sh

sed -r "s/Commit: .*/Commit: $(git log | cut -d ' ' -f 2 | head -n 1)/" DESCRIPTION > out; mv out DESCRIPTION
git add DESCRIPTION
git commit --amend
