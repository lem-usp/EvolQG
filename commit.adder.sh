#!/bin/sh

sed -r "s/Commit: .*/Commit: $(git log | cut -d ' ' -f 2 | head -n 1)/" DESCRIPTION > out; mv out DESCRIPTION
<<<<<<< HEAD
git add DESCRIPTION
git commit --amend
=======
>>>>>>> 8b7535e3ce3b1c157689d40aa995071890422784
